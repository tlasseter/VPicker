package com.Sts.Actions.Wizards.FractureInterpret;

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

public class StsFractureConstructionPanel extends StsJPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private String name = "FractureName";
    private StsColor stsColor;
    private StsWell well;

    private StsModel model = null;

	StsGroupBox panelBox = new StsGroupBox("Define New Fracture");
    StsStringFieldBean fractureName = new StsStringFieldBean();
    StsColorComboBoxFieldBean fractureColor = new StsColorComboBoxFieldBean();
    StsComboBoxFieldBean wellBean = new StsComboBoxFieldBean();

    public StsFractureConstructionPanel(StsWizard wizard, StsWizardStep wizardStep)
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
        fractureName.initialize(this, "name", true, "Name:");
        StsColor[] colors = model.getSpectrum("Basic").getStsColors();
		stsColor = colors[0];
        fractureColor.initializeColors(this, "stsColor", "Color:", colors);
        StsObject[] wells = model.getObjectList(StsWell.class);
        wellBean.initialize(this, "well", "Well", wells);
        gbc.fill = gbc.HORIZONTAL;
		add(panelBox);
        panelBox.gbc.fill = gbc.HORIZONTAL;
		panelBox.add(fractureName);
		panelBox.add(fractureColor);
        panelBox.add(wellBean);
     }


    public String getName() { return name; }
    public StsColor getStsColor() { return stsColor; }
    public void setName(String name)
    {
        this.name = name;
        if(name != "") wizard.enableFinish();
    }
    public void setStsColor(StsColor color) { stsColor = color; }

    public void setWell(Object well) { this.well = (StsWell)well; }
    public Object getWell() { return well; }
}