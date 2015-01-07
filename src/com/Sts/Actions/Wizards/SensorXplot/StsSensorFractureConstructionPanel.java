package com.Sts.Actions.Wizards.SensorXplot;

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

public class StsSensorFractureConstructionPanel extends StsJPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private String name = "FractureName";
    private StsColor stsColor;

    private StsModel model = null;

	StsGroupBox panelBox = new StsGroupBox("Define New Fracture");
    StsStringFieldBean fractureName = new StsStringFieldBean();
    StsColorComboBoxFieldBean fractureColor = new StsColorComboBoxFieldBean();

    public StsSensorFractureConstructionPanel(StsWizard wizard, StsWizardStep wizardStep)
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

//        setLayout(gridBagLayout3);
        model = wizard.getModel();
        fractureName.initialize(this, "name", true, "Name:");
        StsColor[] colors = model.getSpectrum("Basic").getStsColors();
		stsColor = colors[0];
        fractureColor.initializeColors(this, "stsColor", "Color:", colors);

        gbc.fill = gbc.HORIZONTAL;
		add(panelBox);
        panelBox.gbc.fill = gbc.HORIZONTAL;
		panelBox.add(fractureName);
		panelBox.add(fractureColor);
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