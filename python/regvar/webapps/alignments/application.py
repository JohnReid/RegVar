import bx.align.maf
import cStringIO
import os

from flask import Flask, request
ENVSETTINGSVAR = 'ALIGNMENTS_SETTINGS'
app = Flask(__name__)
app.config.from_object('regvar.webapps.alignments.DevelopmentConfig')
if ENVSETTINGSVAR in os.environ:
    app.config.from_envvar(ENVSETTINGSVAR)


def get_maf(genome, alignment, chrom):
    return os.path.join(
        app.config['UCSC_DIR'], 'goldenPath', genome,
        alignment, 'maf', '{0}.maf.bz2'.format(chrom))


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
    app.run(host='0.0.0.0', port=9083)