from cc3d import CompuCellSetup

from chemicalmovengrowSteppables import CellInitializerSteppable
CompuCellSetup.register_steppable(steppable=CellInitializerSteppable(frequency=1))


from chemicalmovengrowSteppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from chemicalmovengrowSteppables import GrowthSteppable

CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))




from chemicalmovengrowSteppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

#from chemicalmovengrowSteppables import scatteredSteppable

#CompuCellSetup.register_steppable(steppable=scatteredSteppable(frequency=1))


CompuCellSetup.run()
