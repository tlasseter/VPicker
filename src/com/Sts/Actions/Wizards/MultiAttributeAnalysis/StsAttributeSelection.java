package com.Sts.Actions.Wizards.MultiAttributeAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsAttributeSelection extends StsWizardStep
{
    StsAttributeSelectionPanel panel;
    StsHeaderPanel header;
    StsMultiAttributeVector vSet = null;

    public StsAttributeSelection(StsWizard wizard)
	{
        super(wizard);
        panel = new StsAttributeSelectionPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Multiple Attribute Analysis");
        header.setSubtitle("Attribute Selection");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#MultiAttributeAnalysis");
        header.setInfoText(wizardDialog,"(1) Select the desired attributes to be used for vector length, color and direction.\n" +
                                        "(2) Specify the length threshold.\n" +
                                        "(3) Press the Finish>> Button");

    }

    public boolean start()
	{
        panel.initialize();
        wizard.enableFinish();
        disableNext();
        return true;
    }

    public boolean end()
    {
        if(vSet != null)
            return true;
        String name = null;
        StsSeismicBoundingBox[] volumes = panel.getSelectedVolumes();
        for(int i=0; i<volumes.length; i++)
        {
        	if(volumes[i] == null)
        	{
        		new StsMessage(wizard.frame, StsMessage.WARNING, "Unable to create multi-vector plot.\n" + 
        				"Must select all three volumes, you can use the same volume in two selections.");
        		return false;
        	}
        }
        //
        // Verify that the volumes have the same bounding box
        //
        /*
        if (volumes.length >= 2)
        {
            if (!volumes[0].sameAs(volumes[1]))
            {
                new StsMessage(wizard.frame, StsMessage.ERROR,
                               "Failed to construct multi-attribute volume. Volumes must have identical bounding box.\n");
                return false;
            }
        }
        */
        //
        // Construct a name if one is not supplied
        //
        name = ((StsMultiAttrAnalysisWizard)wizard).getVolumeName();
        if(name.equals("MAVName") || (name == null) || (name.length() <= 0))
        {
            name = makeVolumeName(volumes);
        }

        // Add the vector set to the object
        vSet = new StsMultiAttributeVector(name, volumes, panel.getNormalizeAzimuth(), panel.getLengthThreshold());
        vSet.addToModel();
        return true;
    }

    public String makeVolumeName(StsSeismicBoundingBox[] volumes)
    {
        String name = "MAV_";
        for(int i=0; i<volumes.length; i++)
        {
        	if(volumes[i] == null)
        		continue;
            if(volumes[i].getName().length() > 8)
                name = name + volumes[i].getName().substring(0,8) + "_";
            else
                name = name + volumes[i].getName() + "_";
        }
        return name;
    }

}
