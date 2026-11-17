
from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)

from pathlib import Path
import json,hashlib,shutil
P=GB_WORK_ROOT
W=P/'wgbs_package'
run=W/'runs/2026-10-14_RNA_DMR_main_supp';out=run/'68_GO_refresh';out.mkdir(exist_ok=True)
src=P/'wgbs/RNA-seq_integration/All_DMR_q05_meth5_merged_overlap_DEGs/DEG_hit_gene_lists_q05_meth5'
for tp,removed,n in [('6hr',{'Pgm1','Sap18b'},1707),('24hr',{'Dynlt1f','Gm5424','Pgm2'},1242)]:
 p=src/f'{tp}_Up_DEGs_hit_by_DMR.txt';v=p.read_text().splitlines();assert len(v)==len(set(v)) and removed<=set(v)
 new=[s for s in v if s not in removed];assert len(new)==n
 (out/f'{tp}_Up_original.txt').write_text('\n'.join(v)+'\n');(out/f'{tp}_Up_corrected.txt').write_text('\n'.join(new)+'\n')
print('Prepared GO gene lists.')
