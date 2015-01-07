package com.Sts.Actions.Wizards.ActiveWell;

import javax.swing.JLabel;
import javax.swing.JScrollPane;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.ProximityAnalysis.StsProximityAnalysisWizard;
import com.Sts.Actions.Wizards.ProximityAnalysis.StsSelectSensors;
import com.Sts.DBTypes.*;
import com.Sts.MVC.StsModel;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectWellsPanel extends StsJPanel
{
    private StsActiveWellWizard wizard;
    private StsSelectWells wizardStep;

    private StsModel model = null;

    transient StsListFieldBean wellListBean;
	transient JScrollPane wellScrollPane = new JScrollPane();
    transient StsButton clearButton = new StsButton();

    public StsSelectWellsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsActiveWellWizard)wizard;
    	this.wizardStep = (StsSelectWells)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
        StsObject[] wells = model.getObjectList(StsWell.class);

        clearButton = new StsButton("Clear Highlights", "Press to clear all the active wells in view.", this, "clearAll");
        wellListBean = new StsListFieldBean(wizard, "selectedWells", null, wells);
        wellScrollPane.getViewport().add(wellListBean, null);
        gbc.fill = gbc.BOTH;
        gbc.gridwidth = 2;
        gbc.weighty = 1.0;
		addEndRow(wellScrollPane);

        gbc.fill = gbc.NONE;
        addEndRow(clearButton);

		wizard.rebuild();
    }

    public void initialize()
    {
    	;
    }

    public void clearAll()
    {
        StsObject[] wells = model.getObjectList(StsWell.class);
        for(int i=0; i<wells.length; i++)
            ((StsWell)wells[i]).setHighlighted(false);
    }

}