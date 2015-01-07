package com.Sts.Actions.Wizards.Fracture;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.StsParameters;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectFractureSetPanel extends StsJPanel
{
    private StsFractureWizard wizard;
    private StsSelectFractureSet wizardStep;

    private StsModel model = null;

	StsGroupBox groupBox = new StsGroupBox("Select or create fracture set");
	StsComboBoxFieldBean fractureSetList;
	StsButton newFractureSetButton;
    StsComboBoxFieldBean timeDepthBean;
    private String zDomainString;

    public StsSelectFractureSetPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsFractureWizard)wizard;
        model = wizard.getModel();
    	this.wizardStep = (StsSelectFractureSet)wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
        fractureSetList = new StsComboBoxFieldBean(wizard, "selectedFractureSet", "Existing Fracture Set", "fractureSetList");
        newFractureSetButton = new StsButton("Define new fracture set", "Define name and color of new fracture set.", wizard, "createNewFractureSet");
        timeDepthBean = new StsComboBoxFieldBean(this, "zDomainString", StsParameters.TD_STRINGS);
        groupBox.gbc.fill = groupBox.gbc.NONE;
        groupBox.gbc.anchor = groupBox.gbc.SOUTH;
        groupBox.add(fractureSetList);
        groupBox.gbc.gridwidth = 2;
        groupBox.add(timeDepthBean);
        groupBox.add(newFractureSetButton);
        gbc.fill = gbc.HORIZONTAL;
		add(groupBox);
    }

    public void initialize()
    {
        Object[] fractureSets = wizard.getFractureSetList();
	    fractureSetList.setListItems(fractureSets);
    }

    public String getZDomainString()
    {
        return model.getProject().getZDomainString();
    }

    public void setZDomainString(String zDomainString)
    {
        if(this.zDomainString == zDomainString) return;
        this.zDomainString = zDomainString;
        model.getProject().setZDomainString(zDomainString);
        fractureSetList.setValueFromPanelObject(wizard);
    }
}
