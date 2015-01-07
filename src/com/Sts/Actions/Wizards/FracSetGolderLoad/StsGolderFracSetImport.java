

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FracSetGolderLoad;

import com.Sts.DBTypes.StsGolderFractureSet;
import com.Sts.IO.StsFile;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.Types.StsGolderFracture;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsStringUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

public class StsGolderFracSetImport
{
	StsModel model;
	StsProject project;

	BufferedReader bufRdr = null;
	String[] propertyNames;
	String[] setNames;
	public String name;
	int currentSetNumber = -1;
	StsGolderFractureSet currentFractureSet;

	static boolean originSet;
	static double xOrigin, yOrigin;

	/** these are the collection of all fractureSets created by all files loaded.
	 *  On completion, these will be added to the model; or not, if aborted during load process */
	ArrayList<StsGolderFractureSet> fractureSets = new ArrayList<StsGolderFractureSet>();

	static public final String BEGIN_STRING = "BEGIN", END_STRING = "END";

	static public final String FORMAT_STRING = "format", PROPERTIES_STRING = "properties", SETS_STRING  = "sets",
			FRACTURES_STRING  = "fracture", TESS_FRACTURES_STRING  = "tessfracture";
	static public final String[] TYPE_KEYWORDS = { FORMAT_STRING, PROPERTIES_STRING , SETS_STRING , FRACTURES_STRING , TESS_FRACTURES_STRING  };
	static public final byte FORMAT = 0, PROPERTIES = 1, SETS  = 21, FRACTURES  = 3, TESS_FRACTURES  = 4;

	static public final String ASCII = "Ascii", XAXIS = "XAxis", EAST = "East", SCALE = "Scale", NO_FRACTURES = "No_Fractures",
			NO_TESS_FRACTURES = "No_TessFractures", NO_NODES = "No_Nodes", NO_PROPERTIES = "No_Properties", NO_NODE_PROPERTIES = "No_Node_Properties";
	static public final String[] FORMAT_KEYWORDS = { ASCII, XAXIS, SCALE, NO_FRACTURES, NO_TESS_FRACTURES, NO_NODES, NO_PROPERTIES, NO_NODE_PROPERTIES, EAST };

	public StsGolderFracSetImport(StsModel model, boolean originSet)
	{
		this.model = model;
		this.project = model.getProject();
		initializeOrigin(originSet);
	}

    public boolean loadFile(StsFile file)
    {
		try
		{
			String filename = file.getFilename();
			name = suffixParser(filename);
			readFile(file);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputWarningException(StsGolderFracSetImport.class, "constructor", e);
			return false;
		}
    }

	private void initializeOrigin(boolean originSet)
	{
		this.originSet = originSet;
		if(originSet)
		 {
			 xOrigin = project.getXOrigin();
			 yOrigin = project.getYOrigin();
		 }
	}

	public void addSetsToProjectAndModel(StsModel model)
	{
		StsProject project = model.getProject();
		for(StsGolderFractureSet fractureSet : fractureSets)
		{
			fractureSet.finish();
			project.addToProject(fractureSet, StsProject.TD_DEPTH);
			project.adjustBoundingBoxes(true, false);
			fractureSet.addToModel();
		}
	}

	//TODO use one of the general parser classes when that stuff is merged (realtime branch)
	static String suffixParser(String filename)
	{
		String[] tokens = StsStringUtils.getTokens(filename, ".");
		return tokens[0];
	}

    public boolean readFile(StsFile file)
    {
        try
        {
            bufRdr = new BufferedReader(new FileReader(file.getPathname()));
            while (true)
            {
				String[] tokens = readTokens(bufRdr);
				if(tokens == null) break;
				if(tokens.length == 0) continue;
				if(tokens[0].equalsIgnoreCase(BEGIN_STRING))
				{
					if(tokens[1].equalsIgnoreCase(FORMAT_STRING))
						readFormats(bufRdr);
					else if(tokens[1].equalsIgnoreCase(SETS_STRING))
						readSets(bufRdr);
					else if(tokens[1].equalsIgnoreCase(PROPERTIES_STRING))
						readProperties(bufRdr);
					else if(tokens[1].equalsIgnoreCase(FRACTURES_STRING))
						readFractures(bufRdr);
					else if(tokens[1].equalsIgnoreCase(TESS_FRACTURES_STRING))
						readTessFractures(bufRdr);
				}
			}
			return true;
        }
         catch(Exception e)
        {
            StsMessageFiles.logMessage("Failed to find time column in file: " + file.getFilename());
            return false;
        }
    }

	private void readFormats(BufferedReader bufRdr)
	{
		while(true)
		{
			String[] tokens = readTokens(bufRdr);
			if(tokens == null) return;
			if(tokens[0].equalsIgnoreCase(END_STRING)) return;
		}
	}

	private void readSets(BufferedReader bufRdr)
	{
		ArrayList<String> setNames = new ArrayList<String>();
		while(true)
		{
			String[] tokens = readTokens(bufRdr);
			if(tokens == null) return;
			if(tokens[0].equalsIgnoreCase(END_STRING)) break;
			setNames.add(tokens[0]);
		}
		this.setNames = setNames.toArray(new String[0]);
	}

	private void readProperties(BufferedReader bufRdr)
	{
		ArrayList<String> propertyNames = new ArrayList<String>();
		while(true)
		{
			String[] tokens = readTokens(bufRdr);
			if(tokens == null) return;
			if(tokens[0].equalsIgnoreCase(END_STRING)) break;
			propertyNames.add(tokens[0]);
		}
		this.propertyNames = propertyNames.toArray(new String[0]);
	}

	private void readFractures(BufferedReader bufRdr)
	{
		while(true)
		{
			String[] tokens = readTokens(bufRdr);
			if(tokens == null) return;
			if(tokens[0].equalsIgnoreCase(END_STRING)) return;
			readFracture(tokens);
		}
	}

	private void readFracture(String[] tokens)
	{
		try
		{
			StsGolderFracture fracture = StsGolderFracture.fractureConstructor(tokens);
			int nVertices = Integer.parseInt(tokens[1]);
			currentSetNumber = Integer.parseInt(tokens[2]);
			fracture.setNumber = currentSetNumber;
			float[][] vertices = new float[nVertices][3];
			for(int n = 0; n < nVertices; n++)
			{
				tokens = readTokens(bufRdr);
				if(tokens == null) return;
				vertices[n] = computeVertexXYZ(tokens);
			}
			fracture.addVertices(vertices);
			tokens = readTokens(bufRdr);
			addNormal(fracture, tokens);
			addFractureToSet(fracture);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "readFractures", e);
		}
	}

	private void readTessFractures(BufferedReader bufRdr)
	{
		while(true)
		{
			String[] tokens = readTokens(bufRdr);
			if(tokens == null) return;
			if(tokens[0].equalsIgnoreCase(END_STRING)) return;
			readTessFracture(tokens);    
		}
	}
	
	private void readTessFracture(String[] tokens)
	{
		try
		{
			int nVertices = Integer.parseInt(tokens[1]);
			int nFaces = Integer.parseInt(tokens[2]);
			currentSetNumber = Integer.parseInt(tokens[3]);
			float[][] allVertices = new float[nVertices][3];
			for(int n = 0; n < nVertices; n++)
			{
				tokens = readTokens(bufRdr);
				if(tokens == null) return;
				allVertices[n] = computeVertexXYZ(tokens);
			}
			for(int n = 0; n < nFaces; n++)
			{
				tokens = readTokens(bufRdr);
				StsGolderFracture fracture = StsGolderFracture.tessFractureConstructor(tokens, allVertices);
				fracture.setNumber = currentSetNumber;
				addFractureToSet(fracture);
			}			
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "readFractures", e);
		}
	}

	private void addNormal(StsGolderFracture fracture, String[] tokens)
	{
		Integer.parseInt(tokens[0]); // normal has index of 0
		float[] normal = new float[3];
		normal[0] = Float.parseFloat(tokens[1]);
		normal[1] = Float.parseFloat(tokens[2]);
		normal[2] = -Float.parseFloat(tokens[3]);
		fracture.addNormal(normal);
	}

	float[] computeVertexXYZ(String[] tokens)
	{
		int nVertex = Integer.parseInt(tokens[0]);
		double x = Double.parseDouble(tokens[1]);
		double y = Double.parseDouble(tokens[2]);
		float z = -Float.parseFloat(tokens[3]);
		 if(!originSet)
		 {
			 xOrigin = x;
			 yOrigin = y;
			 x = 0.0;
			 y = 0.0;
			 originSet = true;
			 project.setOrigin(xOrigin, yOrigin);
		 }
		else
		 {
			 x -= xOrigin;
			 y -= yOrigin;
		 }

		float[] xy = model.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(x, y);
		return new float[] { (float)x, (float)y, z };
	}

	void addFractureToSet(StsGolderFracture fracture)
	{
		String fractureSetName = setNames[fracture.setNumber-1];
		if(currentFractureSet == null || !currentFractureSet.name.equals(fractureSetName))
			currentFractureSet = getCreateFractureSet(fractureSetName);
		currentFractureSet.addFracture(fracture);
	}

	public StsGolderFractureSet getCreateFractureSet(String fractureSetName)
	{
		for(StsGolderFractureSet fractureSet : fractureSets)
		{
			if(fractureSet.name.equals(fractureSetName))
			{
				currentFractureSet = fractureSet;
				return fractureSet;
			}
		}
		StsGolderFractureSet fractureSet = StsGolderFractureSet.constructor(fractureSetName, propertyNames, xOrigin, yOrigin);
		fractureSets.add(fractureSet);
		return fractureSet;
	}

	private String getCurrentSetName()
	{
		return setNames[currentSetNumber-1];
	}

    public String[] readTokens(BufferedReader bufRdr)
    {
    	try
    	{
    		String line = bufRdr.readLine();
			if(line == null) return null;
			String[] tokens = StsStringUtils.getTokens(line);
			if(tokens == null)  return new String[0];
    		return tokens;
    	}
    	catch(Exception ex)
    	{
    		return null;
    	}
    }
}
