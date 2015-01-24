select L.id, M.pdbFile from  LigandSuperRelationEJB as L left join MoadletEJB as M on (L.moadlet=M.id)
