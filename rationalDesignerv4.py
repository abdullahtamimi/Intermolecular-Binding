#!/usr/bin/env python3

from pymol import cmd,stored

negAA={"ASP","GLU"};
posAA={"LYS","ARG","HIS"}
ringAA={"TYR","TRP","PHE"}
covalentAA={"CYS"}
hbondAA={"SER","THR","TYR","ASN","GLN","HIS","MET","CYS","TRP"}
nonPolarAA={"GLY","ALA","ILE","LEU","VAL","PHE","PRO","TRP","MET"}

def pairFind(receptor,ligand,states=0,
                   penaltyDist="4.4",ionicDist="3",hDist="3.5",
                   covalentDist="1.8",ringDist="3.5",nonPolarDist='2.5'):

    '''Takes input of two protein selections in PyMOL, a 'receptor' and 'ligand',
    as well as modifiable distances for the different types of bonds this script will look at.
    These distances will not be used in this function, but are included 
    
    Will look at all amino acids within a certain distance cutoff about the area in contact with
    receptor and ligand, iterate through these residues in perspective of the receptor,
    and find amino acids in the ligand that are a set distance away.
    
    Identifiying info from these will be stored in a dictionary, where receptor amino acids,
    resi number, and chain are the identifying keys, and conjugate amino acids are values.
    This dictionary will be passed on to the next function called'''
    
    pairDict=dict()
    stored.receptorContacts=set()
    
    if states==0:
 
	cmd.split_states(ligand)
	cmd.select('contact',"("+receptor+" and ("+ligand+"_* around 4))")
	cmd.delete(ligand+'_*')
	
	cmd.iterate("contact",'stored.receptorContacts.add((resn,"resi "+resi,"chain "+chain))')

	for items in stored.receptorContacts:
			
	    receptorResId=str(items[0]+":"+items[1]+":"+items[2])
			
	    resnPair={'receptorAA':set()}
                            
	    cmd.select('residuePair','byRes '+ligand+' within 4 of '+items[1])
	    cmd.iterate('residuePair','receptorAA.add(resn+":"+"resi "+resi+":"+"chain "+chain)',space=resnPair)
	    cmd.delete("residuePair")
	    pairDict[receptorResId]=resnPair['receptorAA']
	    del resnPair['receptorAA']
	     
        pairClassify(receptor,ligand,pairDict,penaltyDist,
                     ionicDist,hDist,covalentDist,ringDist,nonPolarDist)
      
def pairClassify(receptor,ligand,pairDict,
                 penaltyDist,ionicDist,hDist,
                 covalentDist,ringDist,nonPolarDist):

    '''input is receptor,ligand,pairDict generated in the pairFind function, and
    distance cutoffs for each type of pairing interaction.

    Will iterate through pairDict and determine/predict if amino acid residues in
    the receptor are forming hydrogen/ionic/covalent/nonpolar bonds with residues in the
    ligand based off distance.  Will also check for 'penalty' bonds, which are amino acid
    residues with like charge that are a close proximity to eachother.

    The amino acid residue pairings will be stored in a dictionary called bondDict,
    with bond type as key, and pairings as values, in a set.  This will be passed on to
    two other functions: pyPrinter which creates interactable objects in PyMOL, and txtPrinter,
    which prints to a text document pairing information.'''

    bondDict=dict()
    hbondSet=set()
    ionicSet=set()
    ringSet=set()
    penaltySet=set()
    covalentSet=set()
    nonPolarSet=set()                    
       
    for receptorResidue,ligandResidue in pairDict.iteritems():
        
        receptorList=receptorResidue.split(":");
	receptorName=receptorList[0]
	
	for ligandHeader in ligandResidue:
        
    	    ligandList=ligandHeader.split(":")
	    ligandName=ligandList[0]
	    			
	    if ligandName in nonPolarAA and receptorName in nonPolarAA:
                if cmd.select("test", ligandList[1]+" and "+ligandList[2]+" within "+nonPolarDist\
                              +" of "+receptorList[1]+" and "+receptorList[2]):
		    nonPolarSet.add(receptorResidue+"%"+ligandHeader)
				    
	    if ligandName in negAA and receptorName in negAA:
		if cmd.select("test","elem o and "+ligandList[1]+" and "+ligandList[2]+" within "+penaltyDist\
                              +" of "+"elem o and "+receptorList[1]+" and "+receptorList[2]):
		    penaltySet.add(receptorResidue+"%"+ligandHeader)
					
	    if ligandName in posAA and receptorName in posAA:
		if cmd.select("test","elem n and "+ligandList[1]+" and "+ligandList[2]+" within "+penaltyDist\
                              +" of "+"elem n and "+receptorList[1]+" and "+receptorList[2]):
		    penaltySet.add(receptorResidue+"%"+ligandHeader)
					
	    if ligandName in ringAA and receptorName in ringAA:
		if cmd.select("test",ligandList[1]+" and "+ligandList[2]+" within "+ringDist\
                              +" of "+receptorList[1]+" and "+receptorList[2]):
		    ringSet.add(receptorResidue+"%"+ligandHeader)
		    		
	    if ligandName in negAA and receptorName in posAA:
		if cmd.select("test","elem o and "+ligandList[1]+" and "+ligandList[2]+" within "+ionicDist\
                              +" of "+"elem n and "+receptorList[1]+" and "+receptorList[2]):
		    ionicSet.add(receptorResidue+"%"+ligandHeader)
					  
            if ligandName in posAA and receptorName in negAA:
		if cmd.select("test","elem n and "+ligandList[1]+" and "+ligandList[2]+" within "+ionicDist\
                              +" of "+"elem o and "+receptorList[1]+" and "+receptorList[2]):
		    ionicSet.add(receptorResidue+"%"+ligandHeader)

            if ligandName in hbondAA and receptorName in hbondAA:
		if cmd.select("test", ligandList[1]+" and "+ligandList[2]+" within "+hDist\
                              +" of "+receptorList[1]+" and "+receptorList[2]):
		    hbondSet.add(receptorResidue+"%"+ligandHeader)
							
            if ligandName in covalentAA and receptorName in covalentAA:
		if cmd.select("test","elem s and "+ligandList[1]+" and "+ligandList[2]+" within "+covalentDist\
                              +" of "+"elem s and "+
                              receptorList[1]+" and "+receptorList[2]):
		    covalentSet.add(receptorResidue+"%"+ligandHeader)
		    
    bondDict['ionic']=ionicSet
    bondDict["hbond"]=hbondSet
    bondDict["aromatic"]=ringSet
    bondDict["penalty"]=penaltySet
    bondDict['covalent']=covalentSet
    bondDict['nonpolar']=nonPolarSet
    cmd.delete('test')

    
    
    pyPrinter(receptor,ligand,bondDict)
    txtPrinter(receptor,ligand,bondDict)
             
def pyPrinter(receptor,ligand,bondDict):

    '''Input to this function is receptor, ligand and bondDict
    formed from pairClassify.

    Will form new objects in PyMOL for each type of bonding interaction
    detected. Does this by iterating through bondDict and concatenating
    a string to be used in PyMOL to make each object selection in PyMOL'''

    hBondPymol=''
    ionicPymol=''
    aromaticPymol=''
    penaltyPymol=''
    covalentPymol=''
    nonPolarPymol=''

    pyObjectDict={'hbond':hBondPymol,'ionic':ionicPymol,'aromatic':aromaticPymol,'penalty':penaltyPymol,
                 'covalent':covalentPymol,'nonpolar':nonPolarPymol}

    for bondType,members in bondDict.iteritems():
        if members:
            for eachBond in members:
                splitBond=eachBond.split("%")
		receptorRes=splitBond[0].split(":")
		ligandRes=splitBond[1].split(":")
		
		for individualResidue in splitBond:
                    parseList=individualResidue.split(":")

                    for bondName in pyObjectDict:
                        if bondType==bondName:
                            #print(parseList[1])
                            pyObjectDict[bondName]+="(" + parseList[1] + " and " + parseList[2] + ")" + " + "         

    cmd.delete("hbond"),cmd.delete("ionic"),
    cmd.delete("penalty");cmd.delete("aromatic"),
    cmd.delete("covalent");cmd.delete("nonpolar")
    
    for objName in pyObjectDict:
        item=pyObjectDict[objName]
        if item.endswith(' + '):
            item=item[:-3]
            cmd.select(objName,item)

def txtPrinter(receptor,ligand,bondDict):
    
    '''Input to this function is receptor, ligand and bondDict
    formed from pairClassify

    Opens a file named residue_pairings where it will write
    amino acid pairings from receptor and ligand, including info
    about residue name, number, and chain it belongs to, as well
    as its conjugate amino acid.

    Pairings will be separated by bond type'''
       
    f=open('residue_pairings.txt','w')
		
    print >>f,("{} and {} residue pairings \n").format(receptor,ligand)
	
    for bondType,members in bondDict.iteritems():
        
        if members:

            print >>f,("{}").format(bondType)

            for eachBond in members:
                
                splitBond=eachBond.split("%")
		receptorRes=splitBond[0].split(":")
		ligandRes=splitBond[1].split(":")
		print(receptorRes[1])

		print >>f, ("\t{:12s}{:10s}  {:3s} ----- {:4s}  {:12s}{:10s}").format(receptorRes[1],receptorRes[2],receptorRes[0],
                                                                                      ligandRes[0],ligandRes[1],ligandRes[2])	
    f.close()
      
cmd.extend("pairFind",pairFind)
