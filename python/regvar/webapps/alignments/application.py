import bx.align.maf
import cStringIO
import os
import io
import ete2
import logging

logger = logging.getLogger(__name__)

from flask import Flask, request, send_file
ENVSETTINGSVAR = 'ALIGNMENTS_SETTINGS'
app = Flask(__name__)
#
# Change this to ProductionConfig when finished
app.config.from_object('regvar.webapps.alignments.DevelopmentConfig')
if ENVSETTINGSVAR in os.environ:
    app.config.from_envvar(ENVSETTINGSVAR)


def get_alignment_dir(genome, alignment):
    return os.path.join(
        app.config['UCSC_DIR'], 'goldenPath', genome, alignment)


def get_maf(genome, alignment, chrom):
    return os.path.join(get_alignment_dir(genome, alignment),
                        'maf', '{0}.maf.bz2'.format(chrom))


def get_treefile(genome, alignment, treename):
    return os.path.join(get_alignment_dir(genome, alignment),
                        '{0}.nh'.format(treename))


@app.route('/newick/<genome>/<alignment>/<treename>')
def newick(genome, alignment, treename):
    return open(get_treefile(genome, alignment, treename)).read()


@app.route('/treeview/<genome>/<alignment>/<treename>')
def treeview(genome, alignment, treename):
    imagetype = request.headers.get('ImageType', 'png')
    import tempfile
    # Construct the tree
    tree = ete2.Tree(open(get_treefile(genome, alignment, treename)).read())
    # Choose rendering style
    ts = ete2.TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = False
    # ts.scale = 120  # 120 pixels per branch length unit
    # Choose a temporary file to store the rendered tree in
    treefilename = tempfile.mktemp(suffix='.' + imagetype, prefix='treeview')
    # Render the tree
    tree.render(treefilename, tree_style=ts)
    # Serve the image
    return send_file(io.BytesIO(open(treefilename).read()),
                     mimetype='image/png')


@app.route('/alignment/<genome>/<alignment>/<chrom>/<int:start>/<int:end>')
def alignment(genome, alignment, chrom, start, end):
    chop = bool(int(request.headers.get('Chop', 0)))
    mincols = int(request.headers.get('Min-Cols', 0))
    src = '{0}.{1}'.format(genome, chrom)
    # app.logger.info(request.headers.get('User-Agent'))
    # app.logger.info(chop)
    index = bx.align.maf.MultiIndexed(
        [get_maf(genome, alignment, chrom)],
        keep_open=True,
        parse_e_rows=True,
        use_cache=True)
    # Write MAF into string
    output = cStringIO.StringIO()
    out = bx.align.maf.Writer(output)
    strand = None
    # Find overlap with reference component
    blocks = index.get(src, start, end)
    # Write each intersecting block
    if chop:
        for block in blocks:
            ref = block.get_component_by_src(src)
            slice_start = max(start, ref.get_forward_strand_start())
            slice_end = min(end, ref.get_forward_strand_end())
            sliced = block.slice_by_component(ref, slice_start, slice_end)
            # If the block is shorter than the minimum allowed size, stop
            if mincols and (sliced.text_size < mincols):
                continue
            # If the reference component is empty, don't write the block
            if sliced.get_component_by_src(src).size < 1:
                continue
            # Keep only components that are not empty
            sliced.components = [c for c in sliced.components if c.size > 0]
            # Reverse complement if needed
            if (strand is not None) and (ref.strand != strand):
                sliced = sliced.reverse_complement()
            # Write the block
            out.write(sliced)
    else:
        for block in blocks:
            out.write(block)
    result = output.getvalue()
    output.close()
    # Close output MAF
    index.close()
    out.close()
    return result

if __name__ == '__main__':
    app.run(host=app.config['HOST'], port=app.config['PORT'])
