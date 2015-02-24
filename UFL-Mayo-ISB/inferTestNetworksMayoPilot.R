require(synapseClient)
synapseLogin()

res <- synQuery('select name, id from file where projectId=="syn2580853" and center=="UFL-Mayo-ISB" and disease=="Alzheimers Disease" and platform=="IlluminaHiSeq2000" and other=="geneCounts" and other=="Normalized"')

synObj <- synGet(res$file.id)
