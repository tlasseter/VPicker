package com.Sts.Actions.Wizards.GriddedSensorAtts;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsSensorVVolume extends StsWizardStep
{
    StsSensorVVolumePanel panel;
    StsHeaderPanel header;
	StsSensorVirtualVolume virtualVolume = null;

    public StsSensorVVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsSensorVVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Virtual Volume Definition");
        header.setSubtitle("Defining Sensor Virtual Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Specify the name of the volume to be constructed\n" +
                            "(2) Select the accumulation type (attribute) to be computed in each cell of the volume\n" +
                            "(3) If the selected accumulation type is dependent on another attribute, select it.\n" +
                            "(4) Specify the Z interval for the volume to be computed, default is existing project Z step\n" +
                            "(5) Select the shape of the search probe, cylindrical or spherical.\n" +
                            "(6) Define the search radius.\n" +
                            "(7) If the selected shape is cylindrical, specify the height of the cylinder.\n" +
                           "     ***** A virtual volume will be created and set to current volume. Move the Z slider to view *****\n");
    }

    public boolean start()
	{
        panel.initialize();
        return true;
    }

    public boolean end()
    {
    	return true;
    }

    public StsTimeCurve getAttribute()
    {
        return (StsTimeCurve)panel.getAttribute();
    }
    public byte getAccumType()
    {
        return panel.getAccumType();
    }
    public float getXyOffset()
    {
        return panel.getXyOffset();
    }
    public float getZOffset()
    {
        return panel.getZOffset();
    }
    public float getZStep()
    {
        return panel.getZStep();
    }
    public byte getShapeType()
    {
        return panel.getShapeType();
    }
}
