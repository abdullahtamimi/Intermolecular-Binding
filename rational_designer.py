from pymol import cmd,stored

negAA={"ASP","GLU"};
posAA={"LYS","ARG","HIS"}
ringAA={"TYR","TRP","PHE"}
covalentAA={"CYS"}
hbondAA={"SER","THR","TYR","ASN","GLN","HIS","MET","CYS","TRP"}
nonPolarAA={"GLY","ALA","ILE","LEU","VAL","PHE","PRO","TRP","MET"}

def surfaceContact(receptor,ligand,states=0,penaltyDist="3.5",ionicDist="3.2",hDist="3.6",covalentDist="2.5",ringDist="3.9",nonPolarDist='3.8'):
	
	stored.receptorContacts=set();pairDict={};bondDict={}
	if states==0:#looking at one conformation of the receptor/ligand interaction
		
		#creates a superposition of the ligand states
		cmd.split_states(ligand)# splits and orients multiple models and multimers from the biological unit file into a set of single-state molecular objects.
		
		
		cmd.select('contact',"("+receptor+" and ("+ligand+"_* around 4))")#selects the residues in contact between ligand state and receptor(additional state made as a precaution)
		cmd.delete(ligand+'_*')#deletes the state after contact surface is made
	    
		cmd.iterate("contact",'stored.receptorContacts.add((resn,"resi "+resi,"chain "+chain))')#iterate through the contacts object and extracts residue name, and number, and chain, stores in a list called receptorContacts

		for items in stored.receptorContacts:# for every amino acid in the receptor that makes contact with an amino acid on the ligand
			
			receptorResId=str(items[0]+":"+items[1]+":"+items[2])#extracts identifying factors for each residue(name,number,chain), separates with ":"
			
			resnPair={'receptorAA':set()}#makes a dictionary for pymol iteration
			
			cmd.select('residuePair','byRes '+ligand+' within 4 of '+items[1])#create an object called residue pair, which finds every residue ~4A away from every iterated residue
			cmd.iterate('residuePair','receptorAA.add(resn+":"+"resi "+resi+":"+"chain "+chain)',space=resnPair)#adds residue pair to a dictionary object with some arbitrary key
			cmd.delete("residuePair")
			pairDict[receptorResId]=resnPair['receptorAA']#stores appropriate contact residue in receptor with the pairing determined above, replacing the key for the the residue in contact
			del resnPair['receptorAA']#deletes dictionary object with arbitrary key
		
		hbondSet=set();ionicSet=set();ringSet=set();penaltySet=set();covalentSet=set();nonPolarSet=set()
		
		for receptorResidue,ligandResidue in pairDict.iteritems():#extracts residue name from receptor
			receptorList=receptorResidue.split(":");
			receptorName=receptorList[0]
			
			for ligandHeader in ligandResidue:#for every residue paired with the receptor residue
			   	ligandList=ligandHeader.split(":")
				ligandName=ligandList[0]#extracts residue name
				
				if ligandName in nonPolarAA and receptorName in nonPolarAA:
					if cmd.select("test",ligandList[1]+" and "+ligandList[2]+" within "+nonPolarDist+" of "+receptorList[1]+" and "+receptorList[2]):
						nonPolarSet.add(receptorResidue+"%"+ligandHeader)
				    
				if ligandName in negAA and receptorName in negAA:#if two residues are like charges
					if cmd.select("test","elem o and "+ligandList[1]+" and "+ligandList[2]+" within "+penaltyDist+" of "+"elem o and "+receptorList[1]+" and "+receptorList[2]):#ensures paired residues are a short distance apart
						penaltySet.add(receptorResidue+"%"+ligandHeader)#storing residue pairings as one string separated by a % symbol
					
				if ligandName in posAA and receptorName in posAA:#if two residues are like charges
					if cmd.select("test","elem n and "+ligandList[1]+" and "+ligandList[2]+" within "+penaltyDist+" of "+"elem n and "+receptorList[1]+" and "+receptorList[2]):
						penaltySet.add(receptorResidue+"%"+ligandHeader)
						
				if ligandName in ringAA and receptorName in ringAA:#if two residues are aromatic rings
					if cmd.select("test",ligandList[1]+" and "+ligandList[2]+" within "+ringDist+" of "+receptorList[1]+" and "+receptorList[2]):
						ringSet.add(receptorResidue+"%"+ligandHeader)
				
				if ligandName in negAA and receptorName in posAA:#if two residues are opposite in charge
					if cmd.select("test","elem o and "+ligandList[1]+" and "+ligandList[2]+" within "+ionicDist+" of "+"elem n and "+receptorList[1]+" and "+receptorList[2]):
						ionicSet.add(receptorResidue+"%"+ligandHeader)
					  
		        if ligandName in posAA and receptorName in negAA:#if two residues are opposite in charge
					if cmd.select("test","elem n and "+ligandList[1]+" and "+ligandList[2]+" within "+ionicDist+" of "+"elem o and "+receptorList[1]+" and "+receptorList[2]):
						ionicSet.add(receptorResidue+"%"+ligandHeader)

		        if ligandName in hbondAA and receptorName in hbondAA:#if two residues can participate in hydrogen bonding
					if cmd.select("test",ligandList[1]+" and "+ligandList[2]+" within "+hDist+" of "+receptorList[1]+" and "+receptorList[2]):
						hbondSet.add(receptorResidue+"%"+ligandHeader)
							
		        if ligandName in covalentAA and receptorName in covalentAA:#if two cysteines
					if cmd.select("test","elem s and "+ligandList[1]+" and "+ligandList[2]+" within "+"elem s and "+covalentDist+" of "+receptorList[1]+" and "+receptorList[2]):
						covalentSet.add(receptorResidue+"%"+ligandHeader)
				

				
		bondDict['ionic']=ionicSet#stores residue pairings to a dictionary with type of bond as the key
		bondDict["hBond"]=hbondSet
		bondDict["aromatic"]=ringSet
		bondDict["penalty"]=penaltySet
		bondDict['covalent']=covalentSet
		bondDict['nonpolar']=nonPolarSet
		cmd.delete('test')
				
		hBondPymol='';ionicPymol='';aromaticPymol='';penaltyPymol='';covalentPymol='';nonPolarPymol=''#bunch of empty strings to be used to create pyMOL objects
		
		f=open('residue_pairings.txt','w')#opens a file to write pairings to
		
		print >>f,("{} and {} residue pairings \n").format(receptor,ligand)
        
        for bondType,members in bondDict.iteritems():#using the bond pairing dictionary I created, I will parse through each key for type of bond
			
			if members:#if there are pairings in this bond type
				print >>f,("{}").format(bondType)#print the bond type name
			
			for eachBond in members:#look at every individual residue pairing
				
				splitBond=eachBond.split("%")#splits bonds by individual residues involved in each
				receptorRes=splitBond[0].split(":");ligandRes=splitBond[1].split(":")#splits identifying factors
				
				print >>f, ("\t{:12s}{:10s}  {:3s} ----- {:4s}  {:12s}{:10s}").format(receptorRes[1],receptorRes[2],receptorRes[0],ligandRes[0],ligandRes[1],ligandRes[2])#prints in a certain format
				
				for individualResidue in splitBond:#look at every individual residue involved in the pairing
					
					parseList=individualResidue.split(":")#and parse through the format I put them in
					
					if bondType is 'hBond':#for each bond type, store the identifying factors in a way pyMOL can read them
						hBondPymol+="(" + parseList[1] + " and " + parseList[2] + ")" + " + "
					if bondType is 'ionic':
						ionicPymol+="(" + parseList[1] + " and " + parseList[2] + ")" + " + "
					if bondType is 'aromatic':
						aromaticPymol+="(" + parseList[1] + " and " + parseList[2] + ")" + " + "
					if bondType is 'penalty':
						penaltyPymol+="(" + parseList[1] + " and " + parseList[2] + ")" + " + "
					if bondType is 'covalent':
						covalentPymol+="(" + parseList[1] + " and " + parseList[2] + ")" + " + "
					if bondType is 'nonpolar':
						nonPolarPymol+="(" + parseList[1] + " and " + parseList[2] + ")" + " + "
       
        f.close()#closes file
        
        cmd.delete("hbond");cmd.delete("ionic");cmd.delete("penalty");cmd.delete("aromatic");cmd.delete("covalent");cmd.delete("nonpolar")#deletes existing objects			
       
        if hBondPymol.endswith(' + '):#remove the extra " + " sign at the end of each string	
            hBondPymol=hBondPymol[:-3]
            cmd.select("hbond",hBondPymol)#create a pyMOL object for each type of bond
        if ionicPymol.endswith(' + '):
            ionicPymol=ionicPymol[:-3]
            cmd.select("ionic",ionicPymol)
        if penaltyPymol.endswith(' + '):
            penaltyPymol=penaltyPymol[:-3]
            cmd.select("penalty",penaltyPymol)
        if aromaticPymol.endswith(' + '):
            aromaticPymol=aromaticPymol[:-3]
            cmd.select("aromatic",aromaticPymol)
        if covalentPymol.endswith(' + '):
            covalentPymol=covalentPymol[:-3]
            cmd.select("covalent",covalentPymol)
        if nonPolarPymol.endswith(' + '):
            nonPolarPymol=nonPolarPymol[:-3]
            cmd.select("nonpolar",covalentPymol)
             
                

cmd.extend("surfaceContact",surfaceContact)#this is so you can call the surfaceContact in pyMOL by typing in surfaceContact in the command window
