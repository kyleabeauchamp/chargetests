import gaff2xml
import openeye.oechem, openeye.oequacpac, openeye.oeomega

iupac = "ethanol"
#iupac = "ethyl_acrylate"
iupac_with_spaces = iupac.replace("_", " ")

for hydrogens in [True, False]:
    for rms in [-1.0, 0.0, 1E-3, 1.0, 1E3]:
        mol = gaff2xml.openeye.iupac_to_oemol(iupac_with_spaces)
        omega = openeye.oeomega.OEOmega()
        omega.SetIncludeInput(True)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(hydrogens)
        eWindow = 15.0
        omega.SetEnergyWindow(eWindow)
        omega.SetMaxConfs(800)

        if rms >= 0.0:
            omega.SetRMSThreshold(rms)

        omega(mol)

        openeye.oequacpac.OEAssignPartialCharges(mol, openeye.oequacpac.OECharges_AM1BCCSym)
        conf = mol.GetConf(openeye.oechem.OEHasConfIdx(0))

        ofs = openeye.oechem.oemolostream()
        ofs.open("./%s/%s_bayly_%s_%s.mol2" % (iupac, iupac, hydrogens, rms))
        openeye.oechem.OEWriteMolecule(ofs, conf)

mol = gaff2xml.openeye.iupac_to_oemol(iupac_with_spaces)
mol1 = gaff2xml.openeye.get_charges(mol)
conf = mol1.GetConf(openeye.oechem.OEHasConfIdx(0))
ofs = openeye.oechem.oemolostream()
ofs.open("./%s/%s_gaff2xml.mol2" % (iupac, iupac))
openeye.oechem.OEWriteMolecule(ofs, conf)

# Analyze output with the following commands:

# grep "C.3       1 <0>        -0" ethanol/*.mol2

# grep "O.2       1 <0>        -0" ethyl_acrylate/*.mol2

