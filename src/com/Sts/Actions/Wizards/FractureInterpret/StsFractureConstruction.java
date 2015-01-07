package com.Sts.Actions.Wizards.FractureInterpret;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsFractureConstruction extends StsWizardStep
{
    StsFractureConstructionPanel panel;
    StsHeaderPanel header;
    StsTriangulatedFracture[] fractures;
    StsFractureInterpretWizard wizard;

    public StsFractureConstruction(StsFractureInterpretWizard wizard)
    {
        super(wizard);
        this.wizard = wizard;
        panel = new StsFractureConstructionPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 600));
        header.setTitle("Define fracture from sensor events");
        header.setSubtitle("");
        header.setInfoText(wizardDialog,"(1) Specify the desired fracture name.\n" +
                                  "(2) Select the fracture color (may change from object panel after creation).\n" +
                                  "(3) Press the Next>> Button to finish or <<Back to define another fracture.");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/FractureInterpret.html");
    }

    public boolean start()
    {
        return constructFracture();
    }

    public boolean end()
    {
        return true;
    }

    private boolean constructFracture()
    {
    	double[][] xyzPoints = null;
    	StsLstSqPointsPlane plane;
        try
        {
            String name = panel.getName();
            if(name.length() == 0)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR, "Please enter a name for this fracture.");
                return false;
            }

            Object[] sensorObjects = wizard.getSelectedSensors();
            int[] clusterNums = wizard.getUniqueClusterNums();
			int nClusterNums = 1;
			if(clusterNums != null) nClusterNums = clusterNums.length;
            fractures = new StsTriangulatedFracture[nClusterNums];
            //
            // Cluster attribute selected, build a fracture for each cluster.
            //
            if(clusterNums != null)
            {
            	for(int i=0; i<clusterNums.length; i++)
            	{
            		for(int j=0; j<sensorObjects.length; j++)
            			xyzPoints = (double[][])StsMath.arrayAddArray(xyzPoints, ((StsDynamicSensor)sensorObjects[j]).getClusteredXYZPoints(clusterNums[i]));
            	
            		if(xyzPoints == null) continue;
            		final StsColor clusterColor = StsColor.colors32[clusterNums[i%32]];
            		fractures[i] = StsTriangulatedFracture.constructor(xyzPoints, name, clusterColor);
            	}
            }
            //
            // No clusters selected, use all points in selected sensors
            //
            else
            {
        		for(int j=0; j<sensorObjects.length; j++)
        			xyzPoints = (double[][])StsMath.arrayAddArray(xyzPoints, ((StsDynamicSensor)sensorObjects[j]).getXYZDoublePoints());
        	
        		if(xyzPoints == null)  return false;
        		final StsColor clusterColor = StsColor.colors32[0];
        		fractures[0] = StsTriangulatedFracture.constructor(xyzPoints, name, clusterColor);
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructFracture", e);
            return false;
        }
    }

	public boolean finish()
	{
		if(fractures == null) return true;
		for(StsTriangulatedFracture fracture : fractures)
			fracture.addToModel();
		return true;
	}
}