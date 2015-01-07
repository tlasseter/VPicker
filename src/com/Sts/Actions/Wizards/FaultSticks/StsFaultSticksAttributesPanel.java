package com.Sts.Actions.Wizards.FaultSticks;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsFaultSticksAttributesPanel extends StsJPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;
    private StsHorpick horpick;

    private String name = "FaultStickSetName";
    private StsColor stsColor;
    private float scalar = 1.0f;
    private boolean overrideColor = false;

    private StsModel model = null;

	StsGroupBox panelBox = new StsGroupBox("Define Fault Stick Set");
    StsStringFieldBean setName = new StsStringFieldBean();   
    StsColorComboBoxFieldBean setColor = new StsColorComboBoxFieldBean();

    public StsFaultSticksAttributesPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;

        try
        {
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        setName.initialize(this, "name", true, "Name:");
        StsColor[] colors = model.getSpectrum("Basic").getStsColors();
		stsColor = colors[3];
        setColor.initializeColors(this, "stsColor", "Color:", colors);
        gbc.fill = gbc.HORIZONTAL;
		add(panelBox);
        panelBox.gbc.fill = gbc.HORIZONTAL;
		panelBox.add(setName);
		panelBox.add(setColor);
     }


    public String getName() { return name; }
    public StsColor getStsColor() { return stsColor; }
    public void setName(String name)
    {
        this.name = name;
        if(name != "") wizard.enableFinish();
    }
    public void setStsColor(StsColor color) { stsColor = color; }
}