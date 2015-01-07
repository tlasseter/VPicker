package com.Sts.Actions.Wizards.SensorXplot;

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

public class StsSensorFractureConstruction extends StsWizardStep
{
    StsSensorFractureConstructionPanel panel;
    StsHeaderPanel header;
    StsTriangulatedFracture[] fracture;
    StsSensorXplotWizard wizard;

    public StsSensorFractureConstruction(StsSensorXplotWizard wizard)
    {
        super(wizard);
        this.wizard = wizard;
        panel = new StsSensorFractureConstructionPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 600));
        header.setTitle("Define fracture from sensor events");
        header.setSubtitle("");
        header.setInfoText(wizardDialog,"(1) Specify the desired fracture name.\n" +
                                  "(2) Select the fracture color (may change from object panel after creation).\n" +
                                  "(4) Press the Next>> Button to finish or <<Back to define another fracture.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorXplot");
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
    	double[][] xyzPoints;
    	StsLstSqPointsPlane plane;
        try
        {
            String name = panel.getName();
            if(name.length() == 0)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR, "Please enter a name for this fracture.");
                return false;
            }
            StsDynamicSensor sensor = wizard.getPrimeSensor();
            int[] clusterNums = wizard.getClusterNums();
            fracture = new StsTriangulatedFracture[clusterNums.length];

            for(int i=0; i<clusterNums.length; i++)
            {
                xyzPoints = sensor.getClusteredXYZPoints(clusterNums[i]);
                if(xyzPoints == null) continue;
                plane = new StsLstSqPointsPlane(xyzPoints);
                final StsColor clusterColor = wizard.getClusterStsColor(i);
                fracture[i] = StsTriangulatedFracture.constructor(plane, name, clusterColor);
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructFracture", e);
            return false;
        }
    }
}