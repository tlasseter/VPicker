package com.Sts.Actions.Wizards.OSWell;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.OpenSpirit.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;
/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsOSWellSelectPanel extends StsJPanel implements ActionListener
{
    StsJPanel unitFormatBox = new StsGroupBox("Units and Format");
    private byte vUnits, hUnits;
    StsComboBoxFieldBean hUnitsBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean vUnitsBean = new StsComboBoxFieldBean();

    StsCheckbox projectLimitBtn = new StsCheckbox("Limit to Project Bounds", "Limit wells to those within project bounds.");

    StsJPanel loadBox = new StsGroupBox("Load Selection");
    StsCheckbox loadTdBtn = new StsCheckbox("Time-Depth", "Load time-depth data.");
    StsCheckbox loadRefBtn = new StsCheckbox("Markers", "Load marker data.");
    StsCheckbox loadLogBtn = new StsCheckbox("Logs", "Load log data.");
    StsCheckbox[] loadBtns = new StsCheckbox[] { loadLogBtn,  loadTdBtn, loadRefBtn };

    StsOSWellObjectTransferPanel transferPanel;
    public StsProgressPanel progressPanel = StsProgressPanel.constructor();

    StsOSWell[] selectedWells = null;

    private StsOSWellWizard wizard;
    private StsOSWellSelect wizardStep;

    public StsOSWellSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsOSWellWizard)wizard;
        this.wizardStep = (StsOSWellSelect)wizardStep;

//        StsOSWellDatastore datastore = StsOSWellDatastore.getInstance();
        StsOSWellDatastore datastore = this.wizard.getWellDatastore();
        datastore.setProjectLimit(true);
        transferPanel = new StsOSWellObjectTransferPanel(datastore, null, wizard.getModel().win3d, progressPanel, 400, 200);
        try
        {
            constructBeans();
            constructPanel();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void constructBeans()
    {
		hUnitsBean.initialize(this, "horzUnitsString", "Horizontal Units:  ", StsParameters.DIST_STRINGS);
		vUnitsBean.initialize(this, "vertUnitsString", "Vertical Units:  ", StsParameters.DIST_STRINGS);
        projectLimitBtn.setSelected(true);
        projectLimitBtn.addActionListener(this);
	}

    public void initialize()
    {
        hUnits = wizard.hUnits;
        vUnits = wizard.vUnits;
		hUnitsBean.setSelectedItem(StsParameters.DIST_STRINGS[hUnits]);
		vUnitsBean.setSelectedItem(StsParameters.DIST_STRINGS[vUnits]);
     }

    void constructPanel() throws Exception
    {
    	gbc.fill = gbc.BOTH;
    	addEndRow(transferPanel);

        gbc.anchor = gbc.NORTH;
        gbc.fill = gbc.HORIZONTAL;
        addEndRow(projectLimitBtn);

        unitFormatBox.addToRow(vUnitsBean);
        unitFormatBox.addEndRow(hUnitsBean);
        addEndRow(unitFormatBox);

        loadBox.addToRow(loadLogBtn);
        loadBox.addToRow(loadRefBtn);
        loadBox.addEndRow(loadTdBtn);
        addEndRow(loadBox);

        addEndRow(progressPanel);
    }

    public void actionPerformed(ActionEvent e)
    {
        if(e.getSource() == projectLimitBtn)
        {
            StsOSWellDatastore datastore = (StsOSWellDatastore)transferPanel.getDatastore();
            datastore.setProjectLimit(projectLimitBtn.isSelected());
            datastore.runGetProjectObjects(transferPanel, progressPanel);
        }
    }

    public StsOSWell[] getSelectedWells()
    {
//        return (StsOSWell[])transferPanel.getSelectedObjects();
    	// this generated a ClassCastException
    	// lkw - workaround
    	Object[] selected = transferPanel.getSelectedObjects();
    	StsOSWell[] wells = new StsOSWell[selected.length];
    	for (int i=0; i<selected.length; i++)
    	{
    		wells[i] = (StsOSWell)selected[i];
    	}
    	return wells;
    }

    public String getHorzUnitsString() { return StsParameters.DIST_STRINGS[hUnits]; }
    public String getVertUnitsString() { return StsParameters.DIST_STRINGS[vUnits]; }
    public byte getHorzUnits() { return hUnits; }
    public byte getVertUnits() { return vUnits; }
    public void setHorzUnitsString(String unitString)
    {
        hUnits = StsParameters.getDistanceUnitsFromString(unitString);
    }

    public void setVertUnitsString(String unitString)
    {
        vUnits = StsParameters.getDistanceUnitsFromString(unitString);
    }

}
