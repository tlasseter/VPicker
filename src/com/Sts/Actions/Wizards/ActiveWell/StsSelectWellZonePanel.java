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
import com.Sts.Utilities.StsMath;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectWellZonePanel extends StsJPanel
{
    private StsActiveWellWizard wizard;
    private StsSelectWellZone wizardStep;

    private StsModel model = null;

    transient StsListFieldBean perfListBean;
    transient StsComboBoxFieldBean wellBean;
	transient JScrollPane perfScrollPane = new JScrollPane();

    transient StsButton clearAllButton = new StsButton();
    transient StsButton clearCurrentButton = new StsButton();

    public StsSelectWellZonePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsActiveWellWizard)wizard;
    	this.wizardStep = (StsSelectWellZone)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
        clearCurrentButton = new StsButton("Clear Current Active", "Press to clear all active perfs from current well.", this, "clearCurrent");
        clearAllButton = new StsButton("Clear All Active", "Press to clear all active perfs from all wells.", this, "clearAll");

        StsObject[] wells = wizard.model.getObjectList(StsWell.class);
        wizard.setSelectedWell((StsWell)wells[0]);
        StsPerforationMarker[] perfs = ((StsWell)wells[0]).getPerforationMarkers();
        
        perfListBean = new StsListFieldBean(wizard, "selectedPerforation", null, perfs);
        perfScrollPane.getViewport().add(perfListBean, null);
        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0;
        wellBean = new StsComboBoxFieldBean(wizard, "selectedWell", "Well:", wells);
    	addEndRow(wellBean);
        gbc.gridwidth = 2;
        JLabel label = new JLabel("Available Perforations");
        addEndRow(label);
        gbc.fill = gbc.BOTH;
        gbc.gridwidth = 2;
        gbc.weighty = 1.0;
		addEndRow(perfScrollPane);
        
        gbc.gridwidth = 1;
        gbc.weighty = 0.0;        
        gbc.fill = gbc.NONE;
        addToRow(clearCurrentButton);
        addEndRow(clearAllButton);

		wizard.rebuild();
    }

    public void rebuildPerfList(StsWell well)
    {
    	StsPerforationMarker[] perfs = well.getPerforationMarkers();
    	if(perfs == null)
    	{
    		new StsMessage(wizard.frame,StsMessage.ERROR, "There are no perfs available for the selected well.\n" +
    				"Please select a different well.");
    		return;
    	}
        perfListBean.setListItems(perfs);
        perfListBean.updateUI();
    }

    public void initialize()
    {
    	StsObject[] wells = wizard.model.getObjectList(StsWell.class);
        StsWell[] activeWells = new StsWell[wells.length];
        int cnt = 0;
        for(int i=0; i<wells.length; i++)
        {
            if(((StsWell)wells[i]).getHighlighted())
                activeWells[cnt++] = (StsWell)wells[i];
        }
        StsMath.trimArray(activeWells, cnt);

        wizard.setSelectedWell(activeWells[0]);
        wellBean.setListItems(activeWells);
    }

    public void clearAll()
    {
        StsObject[] wells = model.getObjectList(StsWell.class);
        for(int i=0; i<wells.length; i++)
        {
            StsPerforationMarker[] perfs = ((StsWell)wells[i]).getPerforationMarkers();
            if(perfs != null)
            {
                for(int j=0; j<perfs.length; j++)
                    perfs[j].setHighlighted(false);
            }
        }
        perfListBean.clearSelections();
    }

    public void clearCurrent()
    {
        StsWell well = wizard.getSelectedWell();
        if(well == null) return;
        StsPerforationMarker[] perfs = well.getPerforationMarkers();
        if(perfs != null)
        {
            for(int j=0; j<perfs.length; j++)
                perfs[j].setHighlighted(false);
        }
        perfListBean.clearSelections();
    }
}