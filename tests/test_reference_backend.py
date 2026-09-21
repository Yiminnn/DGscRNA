"""Regression contracts for packaging the scientific reference implementation."""
import ast
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "dgscrna" / "reference" / "backend"
SOURCE = ROOT / "research" / "reference_examples" / "backend"


class SchedulerMetadataNormalizer(ast.NodeTransformer):
    """Only the provenance value changed; mathematical statements must be identical."""

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id == "execution_id":
            return ast.parse("os.environ['SLURM_JOB_ID']", mode="eval").body
        return self.generic_visit(node)


def component(path, name):
    tree = ast.parse(path.read_text())
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.dump(SchedulerMetadataNormalizer().visit(node), include_attributes=False)
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == name
                                               for target in node.targets):
            return ast.dump(node, include_attributes=False)
    raise AssertionError(f"Missing {name} in {path}")


class ReferenceBackendTests(unittest.TestCase):
    def test_historical_dl_numerical_body_is_identical(self):
        for name in ("PARAMS", "build_model", "train_cache"):
            with self.subTest(component=name):
                self.assertEqual(component(SOURCE / "legacy_refine.py", name),
                                 component(BACKEND / "legacy_refine.py", name))

    def test_source_provenance_checksums(self):
        manifest = json.loads((BACKEND / "SOURCE_DERIVATION.json").read_text())
        for name, item in manifest["files"].items():
            with self.subTest(file=name):
                self.assertEqual(hashlib.sha256((BACKEND / name).read_bytes()).hexdigest(),
                                 item["packaged_sha256"])
                if "source" in item:
                    self.assertEqual(hashlib.sha256((ROOT / item["source"]).read_bytes()).hexdigest(),
                                     item["source_sha256"])

    def test_frozen_r_arithmetic_blocks_are_retained(self):
        prepare_source = (SOURCE / "prepare_R.R").read_text()
        prepare = (BACKEND / "prepare_R.R").read_text()
        start = "obj<-NormalizeData(obj,normalization.method='LogNormalize'"
        end = "script<-sub('^--file='"
        source_body = prepare_source[prepare_source.index(start):prepare_source.index(end)]
        package_body = prepare[prepare.index(start):prepare.index(end)]
        self.assertEqual(source_body, package_body)
        score_source = (SOURCE / "score_R.R").read_text()
        score = (BACKEND / "score_R.R").read_text()
        for start, end in [("    if(method=='SNN')", "  m<-list(status='score_complete_DL_pending'")]:
            self.assertEqual(score_source[score_source.index(start):score_source.index(end)],
                             score[score.index(start):score.index(end)])

    def test_no_site_specific_backend_paths(self):
        for path in BACKEND.iterdir():
            if path.suffix in {".py", ".R"}:
                self.assertNotIn("/fs/scratch/", path.read_text(), str(path))
                self.assertNotIn("/users/PCON", path.read_text(), str(path))

    def test_terminal_csv_readers_preserve_string_identifiers(self):
        tree = ast.parse((BACKEND / "terminal.py").read_text())
        readers = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                   and isinstance(node.func, ast.Attribute) and node.func.attr == "read_csv"]
        self.assertEqual(len(readers), 3)
        for reader in readers:
            keywords = {item.arg: item.value for item in reader.keywords}
            self.assertEqual(ast.dump(keywords["dtype"]), ast.dump(ast.Name(id="str", ctx=ast.Load())))
            self.assertIs(keywords["keep_default_na"].value, False)

    def test_r_csv_readers_preserve_ids_and_numeric_coordinates(self):
        rscript = os.environ.get("DGSCRNA_RSCRIPT") or shutil.which("Rscript")
        if not rscript:
            self.skipTest("Set DGSCRNA_RSCRIPT to check actual R CSV reader expressions")
        script = r'''
args<-commandArgs(trailingOnly=TRUE)
backend<-args[[1]]
assignment<-function(file,name) {
  found<-list()
  walk<-function(x) {
    if(is.call(x)) {
      if(identical(x[[1]],as.name('<-')) && identical(x[[2]],as.name(name)))
        found[[length(found)+1L]]<<-x
      for(i in seq_along(x)[-1])if(!identical(x[[i]],quote(expr=)))walk(x[[i]])
    } else if(is.expression(x) || is.list(x))
      for(i in seq_along(x))if(!identical(x[[i]],quote(expr=)))walk(x[[i]])
  }
  walk(parse(file=file));stopifnot(length(found)==1L);found[[1]]
}
work<-tempfile('reference-identifiers-');dir.create(work)
ids<-c('001','NA','1e3','TRUE')
fit<-data.frame(cell_id=ids,batch=rep('001',4),stringsAsFactors=FALSE)
write.csv(fit,file.path(work,'cells_fit.csv'),row.names=FALSE)
write.csv(fit,file.path(work,'cells.csv'),row.names=FALSE)
write.csv(data.frame(gene=ids),file.path(work,'genes.csv'),row.names=FALSE)
write.csv(data.frame(cell_id=ids,cluster=c('01','NA','002','0')),
          file.path(work,'clusters.csv'),row.names=FALSE)
z<-matrix(seq_len(8)/10,4,2,dimnames=list(ids,c('UMAP_1','UMAP_2')))
write.csv(z,file.path(work,'UMAP2.csv'))
inp<-prep<-dest<-work
input_paths<-list(canonical_cells=file.path(work,'cells.csv'))
umap_path<-file.path(work,'UMAP2.csv')
read_cache<-function(path,reader,valid)reader(path)
eval(assignment(file.path(backend,'prepare_R.R'),'cells'))
stopifnot(identical(cells$cell_id,ids),identical(cells$batch,rep('001',4)))
eval(assignment(file.path(backend,'prepare_R.R'),'genes'))
stopifnot(identical(genes,ids))
eval(assignment(file.path(backend,'score_R.R'),'cells'))
stopifnot(identical(cells$cell_id,ids))
eval(assignment(file.path(backend,'score_R.R'),'cf'))
stopifnot(identical(cf$cell_id,ids),identical(cf$cluster,c('01','NA','002','0')))
eval(assignment(file.path(backend,'score_R.R'),'embedding'))
stopifnot(identical(rownames(embedding),ids),is.numeric(embedding),identical(embedding,z))
eval(assignment(file.path(backend,'ptc_prepare.R'),'canonical'))
stopifnot(identical(canonical$cell_id,ids))
unlink(work,recursive=TRUE)
cat('R_IDENTIFIER_ROUNDTRIP_OK\n')
'''
        result = subprocess.run([rscript, "--vanilla", "-e", script, str(BACKEND)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R_IDENTIFIER_ROUNDTRIP_OK", result.stdout)


if __name__ == "__main__":
    unittest.main()
