select LigandSuperRelationEJB.id, MoadletEJB.pdbFile, LigandEJB.name, LigandEJB.mw, BindingDataEJB.affinmicromolar, BinEJB.ecNumber, , LigandEJB.smiles from LigandSuperRelationEJB left join LigandEJB on LigandSuperRelationEJB.ligandId = LigandEJB.id left join MoadletEJB on LigandSuperRelationEJB.moadlet = MoadletEJB.id left join BindingDataEJB on LigandSuperRelationEJB.bindingData = BindingDataEJB.id left join SubbinEJB on MoadletEJB.subbin = SubbinEJB.id left join BinEJB on SubbinEJB.bin = BinEJB.id;