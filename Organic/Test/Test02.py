import sys
from openeye import oechem
from openeye import oedocking

def main(argv=[__name__]):
    itf = oechem.OEInterface(InterfaceData)
    oedocking.OEScoreTypeConfigure(itf, "-score")
    
    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, itf.GetString("-receptor")):
        oechem.OEThrow.Fatal("Unable to read receptor")

    imstr = oechem.oemolistream()
    if not imstr.open(itf.GetString("-in")):
        oechem.OEThrow.Fatal("Unable to open input file of ligands")

    omstr = oechem.oemolostream()
    if not omstr.open(itf.GetString("-out")):
        oechem.OEThrow.Fatal("Unable to open out file for rescored ligands")

    scoreType = oedocking.OEScoreTypeGetValue(itf, "-score")
    score = oedocking.OEScore(scoreType)
    score.Initialize(receptor)

    for ligand in imstr.GetOEMols():
        if itf.GetBool("-optimize"):
            score.SystematicSolidBodyOptimize(ligand)
        score.AnnotatePose(ligand)
        sdtag = score.GetName()
        oedocking.OESetSDScore(ligand, score, sdtag)
        oechem.OESortConfsBySDTag(ligand, sdtag, score.GetHighScoresAreBetter())
        oechem.OEWriteMolecule(omstr, ligand)

    return 0

InterfaceData = """
!PARAMETER -receptor
!ALIAS -rec
!TYPE string
!REQUIRED true
!LEGAL_VALUE *.oeb
!LEGAL_VALUE *.oeb.gz
!BRIEF A receptor file the poses pass to the -in flag will be scored against
!END
!PARAMETER -in
!TYPE string
!REQUIRED true
!BRIEF Input molecule file with poses to rescore
!END
!PARAMETER -out
!TYPE string
!REQUIRED true
!BRIEF Rescored molecules will be written to this file
!END
!PARAMETER -optimize
!ALIAS -opt
!TYPE bool
!DEFAULT false
!BRIEF Optimize molecules before rescoring
!END
"""

if __name__ == "__main__":
    sys.exit(main(sys.argv))
