select L.id, M.pdbFile,G.name  from  LigandSuperRelationEJB as L left join MoadletEJB as M on (L.moadlet=M.id)  left join LigandEJB as G on (L.ligandId=G.id) where G.name <> ""
