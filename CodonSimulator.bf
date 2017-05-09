#include "simulation.config";				
				

/* ________________________________________________

   	Functions needed to generate a random tree
   
  __________________________________________________ */

function TreeMatrix2TreeString (doLengths)
{
	treeString = "";
	d=treeString*(seqCount*25);
	p = 0;
	k = 0;
	m = treeNodes[0][1];
	n = treeNodes[0][0];

	while (m)
	{	
		if (m>p)
		{
			if (p)
			{
				d=treeString*",";
			}
			for (j=p;j<m;j=j+1)
			{
				d=treeString*"(";
			}
		}
		else
		{
			if (m<p)
			{
				for (j=m;j<p;j=j+1)
				{
					d=treeString*")";
				}
			}	
			else
			{
				d=treeString*",";
			}	
		}
		if (n<seqCount)
		{
			d=treeString*Format(n+1,0,0);
		}
		if (doLengths>.5)
		{
			d=treeString*(":"+treeNodes[k][2]);
		}
		k=k+1;
		p=m;
		n=treeNodes[k][0];
		m=treeNodes[k][1];
	}

	for (j=m;j<p;j=j+1)
	{
		d=treeString*")";
	}
	
	d=treeString*0;
	fprintf(stdout,"\n ",treeString, "\n\n");
	return treeString;
}


/* ____________________________________________*/

function GenerateSymmetricTree (seqCount, bL)
{
	parentInfo		= {2*(seqCount+1),3};
	
	curCount = seqCount/2;
	
	for (k=0; k<seqCount; k=k+1)
	{
		parentInfo[k][2] = 1;
	}	
	
	curPos   = 0;
	curPlace = seqCount;
	
	while (curCount)
	{
		for (k=0; k<curCount; k=k+1)
		{
			parentInfo[curPos][0]=curPlace; /* parent */
			parentInfo[curPos][1]=1;

			parentInfo[curPos+1][0]=curPlace; /* parent */
			parentInfo[curPos+1][1]=1;
			
			parentInfo[curPlace][2] = parentInfo[curPos][2]+parentInfo[curPos+1][2]+1;
			
			curPos = curPos+2;
			curPlace = curPlace+1;
		}
		curCount = curCount $ 2;
	}
	
	parentInfo[curPlace-1][0]=-1; /* parent */
	parentInfo[curPlace-1][1]=1;
	
	njm = parentInfo >= seqCount;

	treeNodes 		= {2*(seqCount+1),3};
	cladesInfo	    = {seqCount-1,2};
	
	for (i=Rows(treeNodes)-1; i>=0; i=i-1)
	{
		treeNodes[i][0] = njm[i][0];
		treeNodes[i][1] = njm[i][1];
		treeNodes[i][2] = njm[i][2];
	}

	for (i=Rows(cladesInfo)-1; i>=0; i=i-1)
	{
		cladesInfo[i][0] = njm[i][3];
		cladesInfo[i][1] = njm[i][4];
	}
	
	njm = 0;
		
	if (bL > 0)
	{
		for (k=0; k<Rows(treeNodes); k=k+1)
		{
			aNumber = Random(0,1);
			FindRoot (root, Exp(-_z_/bL)-1+aNumber ,_z_, 0, 100000);			
			treeNodes[k][2] = root;
		}
	}
	
	fprintf (stdout, "\nGenerating Tree String\n");
	return TreeMatrix2TreeString (bL>0);
}

/* ____________________________________________*/

function GenerateLadderTree (seqCount, bL)
{
	parentInfo		= {2*(seqCount+1),3};
	
	curCount = seqCount/2;
	
	for (k=0; k<seqCount; k=k+1)
	{
		parentInfo[k][2] = 1;
	}	
	
	
	parentInfo[0][0]=seqCount;
	parentInfo[0][1]=1;
	parentInfo[1][0]=seqCount;
	parentInfo[1][1]=1;
	parentInfo[seqCount][2]=3;

	curPos   = 2;
	curPlace = seqCount+1;

	for (curPos=2; curPos<seqCount; curPos=curPos+1)
	{
		parentInfo[curPos][0]=curPlace; /* parent */
		parentInfo[curPos][1]=1;

		parentInfo[curPlace-1][0] = curPlace;
		parentInfo[curPlace-1][1] = 1;
		parentInfo[curPlace][2] = parentInfo[curPlace-1][2]+2;
		curPlace = curPlace+1;
	}
	
	parentInfo[curPlace-1][0]=-1; /* parent */
	parentInfo[curPlace-1][1]=1;
	
	njm = parentInfo >= seqCount;

	treeNodes 		= {2*(seqCount+1),3};
	cladesInfo	    = {seqCount-1,2};
	
	for (i=Rows(treeNodes)-1; i>=0; i=i-1)
	{
		treeNodes[i][0] = njm[i][0];
		treeNodes[i][1] = njm[i][1];
		treeNodes[i][2] = njm[i][2];
	}

	for (i=Rows(cladesInfo)-1; i>=0; i=i-1)
	{
		cladesInfo[i][0] = njm[i][3];
		cladesInfo[i][1] = njm[i][4];
	}
	
	njm = 0;
		
	if (bL > 0)
	{
		for (k=0; k<Rows(treeNodes); k=k+1)
		{
			aNumber = Random(0,1);
			FindRoot (root, Exp(-_z_/bL)-1+aNumber ,_z_, 0, 100000);			
			treeNodes[k][2] = root;
		}
	}
	
	fprintf (stdout, "\nGenerating Tree String\n");
	return TreeMatrix2TreeString (bL>0);
}


/* ____________________________________________*/

function GenerateARandomTree (seqCount, bL)
{
	parentInfo		= {2*(seqCount+1),3};
	canUse			= {seqCount, 1};
	columnIndex 	= {seqCount, 1};
	
	for (idx = 0; idx < seqCount; idx = idx + 1)
	{
		parentInfo [idx][2] = 1;
		columnIndex[idx]    = idx;
	}

	cladesMade 	= 1;
	toDo		= seqCount-1;
	
	sThresh = Sqrt(seqCount/2)$1;
	
	switchingThresh = seqCount - sThresh;
	
	echoStep    = toDo$100 + 1;
	lastDone	= 1;

	while (cladesMade < seqCount)
	{
		if (cladesMade < switchingThresh)
		{
			m = Random (0,seqCount-0.00000000001);
			while (canUse[m] < 0)
			{
				m = Random (0,seqCount-0.00000000001);		
			}
			canUse[m] = -1;
			
			n = Random (0,seqCount-0.00000000001);
			while (canUse[n] < 0)
			{
				n = Random (0,seqCount-0.00000000001);		
			}
		}
		else
		{
			if (cladesMade == switchingThresh)
			{
				newLength = 0;
				m = {sThresh+1,1};
				n = {sThresh+1,1};
				for (iCount = 0; iCount < seqCount; iCount = iCount + 1)
				{
					if (canUse[iCount] == 0)
					{
						m[newLength] = columnIndex[iCount];
						newLength = newLength + 1;
					}
				}
				canUse = n;
				columnIndex = m;
				cladesMade2 = 0.000000000001;
			}
			
			m = Random(0, newLength-cladesMade2)$1;
			
			iCount = -1;
			i = 0;
			
			while (iCount < m)
			{
				if (canUse[i] == 0)
				{
					iCount = iCount + 1;
				}
				i = i+1;
			}
			
			n = m;
			
			for (;n==m;)
			{
				n = Random(0, newLength-cladesMade2)$1;
			}
			
			m = i-1;
			
			iCount = -1;
			i = 0;
			
			while (iCount < n)
			{
				if (canUse[i] == 0)
				{
					iCount = iCount + 1;
				}
				i = i+1;
			}
			n = i-1;
			
			canUse[m] = -1;
			
			cladesMade2 = cladesMade2 + 1;
		}
			
		
		minIndex = n;
		
		m = columnIndex[m];
		n = columnIndex[n];
		
		k = seqCount+cladesMade-1;

		parentInfo[n][0]=k; /* parent */
		parentInfo[n][1]=1;

		parentInfo[m][0]=k; /* parent */
		parentInfo[m][1]=1;

		parentInfo[k][2] = parentInfo[n][2]+parentInfo[m][2]+1;
		

		columnIndex  [minIndex]  = k;
		
		cladesMade = cladesMade + 1;
		toDo = toDo-1;
		
		if (cladesMade-lastDone > echoStep)
		{
			/*
			fprintf (stdout, cladesMade, " nodes selected (", 100*cladesMade/(seqCount-1), "%)\n");
			*/
			lastDone = cladesMade;
		}
	}
	
	parentInfo[k][0]=-1; /* parent */
	parentInfo[k][1]=1;
	
	njm = parentInfo >= seqCount;

	treeNodes 		= {2*(seqCount+1),3};
	cladesInfo	    = {seqCount-1,2};
	
	for (i=Rows(treeNodes)-1; i>=0; i=i-1)
	{
		treeNodes[i][0] = njm[i][0];
		treeNodes[i][1] = njm[i][1];
		treeNodes[i][2] = njm[i][2];
	}

	for (i=Rows(cladesInfo)-1; i>=0; i=i-1)
	{
		cladesInfo[i][0] = njm[i][3];
		cladesInfo[i][1] = njm[i][4];
	}
	
	njm = 0;
		
	if (bL > 0)
	{
		for (k=0; k<Rows(treeNodes); k=k+1)
		{
			aNumber = Random(0,1);
			FindRoot (root, Exp(-_z_/bL)-1+aNumber ,_z_, 0, 100000);			
			treeNodes[k][2] = root;
		}
	}
	
	fprintf (stdout, "\nGenerating Tree String\n");
	
	return TreeMatrix2TreeString (bL>0);
}

/* MAIN PORTION BEGINS HERE */

  

/* compute branch conversion factor */



R  = 1;

synRate    = 1;

blf = 0;

for (k=0; k<ModelMatrixDimension; k=k+1)
{
	for (k2=0; k2<ModelMatrixDimension; k2=k2+1)
	{
		if (k2!=k)
		{
			blf = blf +  vectorOfFrequencies[k]*MG94custom[k][k2];
		}
	}
}

brLen = blf/3;

if (treeType)
{
	if (treeType == 1)
	{
		treeString = GenerateARandomTree(seqs,branchLength/brLen);
	}
	else
	{
		if (treeType == 2)
		{
			treeString = GenerateLadderTree(seqs,branchLength/brLen);
		}
	}
}
else
{
	treeString = GenerateSymmetricTree(seqs,branchLength/brLen);
}

treeNodes 			= 0;
MESSAGE_LOGGING 	= -1;
ACCEPT_ROOTED_TREES = 1;

Tree   					T = treeString;

if (treeType == 3)
{
	branchNames = BranchName (T,-1);
	for (it = 0; it < Columns (branchNames); it = it+1)
	{
		ExecuteCommands ("T."+branchNames[it]+".synRate=T."+branchNames[it]+".synRate/brLen;");
	}
}

DATAFILE_TREE 			= Format (T,0,0);

charactersUniversalCode = {{"A","C","G","T"}
			  			   {"3",GeneticCodeExclusions,"",""}};
			  			   
SetDialogPrompt ("Save simulated data to:");
fprintf 		(PROMPT_FOR_FILE, CLEAR_FILE);
fPath 		= LAST_FILE_PATH;

fprintf (stdout, "\n");

for (it = 0; it < iterates; it = it+1)
{		  
	for (k=0; k<Columns (omegas); k=k+1)
	{
		R = omegas[2][k];
		if (k>0)
		{
			ClearConstraints (T2);
			Tree 		   T2 = treeString;
			R2 = omegas[0][k];
			ReplicateConstraint ("this1.?.synRate:=R2__*this2.?.synRate",T2,T);
			DataSet sim   = Simulate 	(T2,vectorOfFrequencies,charactersUniversalCode,omegas[1][k],0);
			DataSet bigDS = Concatenate (bigDS, sim);
		}
		else
		{
			DataSet bigDS = Simulate 	  (T,vectorOfFrequencies,charactersUniversalCode,omegas[1][k],0);
		}
	}			  				

	DataSetFilter 	bigFilter = CreateFilter (bigDS,1);

	IS_TREE_PRESENT_IN_DATA = 1;
	fName = fPath + "." + it;
  

  fprintf       (stdout, vectorOfFrequencies);
  
	fprintf 		(fName, CLEAR_FILE,bigFilter);
	fprintf			(stdout, "Done with iterate ", it+1, "/", iterates, "\n");
}
