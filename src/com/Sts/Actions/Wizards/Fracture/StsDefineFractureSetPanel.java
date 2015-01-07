package com.Sts.Actions.Wizards.Fracture;

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

public class StsDefineFractureSetPanel extends StsJPanel
{
    private StsFractureWizard wizard;
    private StsDefineFractureSet wizardStep;
 
    private StsModel model = null;

    static int newNameIndex = 0;

    StsGroupBox panelBox = new StsGroupBox("Define New Fracture Set");
    StsStringFieldBean nameBean;
    StsColorComboBoxFieldBean colorBean;
	StsGroupBox groupBox = new StsGroupBox("Use same Z locations as...");
	StsComboBoxFieldBean fractureSetList;
	
    public StsDefineFractureSetPanel(StsFractureWizard wizard, StsDefineFractureSet wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
        nameBean = new StsStringFieldBean(wizard, "fractureSetName", true, "Name:");
        fractureSetList = new StsComboBoxFieldBean(wizard, "selectedFractureSet", "Fracture Set:", wizard.getFractureSetListWithNone());        
        StsSpectrum spectrum = wizard.getModel().getSpectrum("Basic");
        colorBean = new StsColorComboBoxFieldBean(wizard, "fractureSetColor", "Color:", spectrum.getStsColors());
        panelBox.gbc.fill = gbc.HORIZONTAL;
        panelBox.gbc.gridwidth = 1;
		panelBox.add(nameBean);
		panelBox.add(colorBean);
		groupBox.add(fractureSetList);
		panelBox.gbc.gridwidth = 3;
		panelBox.add(groupBox);
        gbc.fill = gbc.HORIZONTAL;
        add(panelBox);
    }

    public void initialize()
    {
        this.initializeBeans();
        Object[] fractureSets = wizard.getFractureSetListWithNone();
	    fractureSetList.setListItems(fractureSets);        
    }
}
