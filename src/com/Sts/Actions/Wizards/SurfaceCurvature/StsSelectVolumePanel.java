package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectVolumePanel extends StsJPanel
{
    private StsVolumeCurvatureWizard wizard;
    private StsSelectVolume wizardStep;
    boolean patchVolSelected = false;
    boolean seisVolSelected = true;

    private StsModel model = null;

    StsComboBoxFieldBean seisVolumeBean;
    StsComboBoxFieldBean patchVolumeBean;
    StsRadioButtonFieldBean selectSeisBean;
    StsRadioButtonFieldBean selectPatchBean;

    public StsSelectVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsVolumeCurvatureWizard)wizard;
    	this.wizardStep = (StsSelectVolume)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
        StsObject[] seisVolumeObjects = model.getObjectList(StsSeismicVolume.class);
        StsObject[] patchVolumeObjects = model.getObjectList(StsPatchVolume.class);
        seisVolumeBean = new StsComboBoxFieldBean(wizard, "selectedVolume", "Seismic Volume:", seisVolumeObjects);
        
        StsGroupBox seisVolBox = new StsGroupBox("Select Seismic Volume");
        StsGroupBox patchVolBox = new StsGroupBox("Select Patch Volume");
        seisVolBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        patchVolBox.gbc.fill = GridBagConstraints.HORIZONTAL;

        selectSeisBean = new StsRadioButtonFieldBean();
        selectSeisBean.initialize(this, "selectSeisVolume", "");
        selectSeisBean.setSelected(true);
        selectSeisBean.addActionListener();
        seisVolBox.addToRow(selectSeisBean);
        seisVolBox.addToRow(seisVolumeBean);
         
        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0;
        addEndRow(seisVolBox);
        if (patchVolumeObjects.length > 0)
        {
        	JLabel label = new JLabel("OR");
        	gbc.fill = gbc.CENTER;
        	addEndRow(label);
            patchVolumeBean = new StsComboBoxFieldBean(wizard, "patchVolume", "Patch Volume:", patchVolumeObjects);
            patchVolumeBean.setSelectedItem(patchVolumeObjects[0]);
            selectPatchBean = new StsRadioButtonFieldBean();
            selectPatchBean.initialize(this, "selectPatchVolume", "");
            selectPatchBean.setSelected(false);
            selectPatchBean.addActionListener();
            patchVolBox.addToRow(selectPatchBean);
	        patchVolBox.addToRow(patchVolumeBean);

        	gbc.fill = gbc.HORIZONTAL;
        	addEndRow(patchVolBox);
        }
        
        
//     	addEndRow(seisVolumeBean);
//        
//    	addEndRow(patchVolumeBean);

		wizard.rebuild();
    }
    
    public boolean getSelectSeisVolume()
    {
    	//patchVolSelected = false;
    	//seisVolSelected = true;
    	return true;
    }
    
    public void setSelectSeisVolume(boolean isSelectd)
    {
    	if (isSelectd)
    	{
    		patchVolSelected = false;
    		seisVolSelected = true;
    		selectPatchBean.setSelected(false);
    	}
    }
    
    public boolean getSelectPatchVolume()
    {
    	//patchVolSelected = true;
    	//seisVolSelected = false;
    	return true;
    }
    
    public void setSelectPatchVolume(boolean isSelectd)
    {
    	if (isSelectd)
    	{
    		patchVolSelected = true;
    		seisVolSelected = false;
    		selectSeisBean.setSelected(false);
    	}
    }
}