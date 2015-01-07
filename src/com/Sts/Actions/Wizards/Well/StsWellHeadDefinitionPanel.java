package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.IO.StsAbstractFile;
import com.Sts.UI.Beans.*;

import javax.swing.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsWellHeadDefinitionPanel extends StsJPanel {
    private StsWellWizard wizard;
    private StsWellHeadDefinition wizardStep;

    JPanel jPanel1 = new JPanel();
    StsGroupBox typeBox = new StsGroupBox("Specify the Wellhead Coordinates");
    StsDoubleFieldBean xBean = new StsDoubleFieldBean();
    StsDoubleFieldBean yBean = new StsDoubleFieldBean();
    StsDoubleFieldBean zBean = new StsDoubleFieldBean();
    StsBooleanFieldBean kbApplyBean = new StsBooleanFieldBean();

    StsGroupBox limitBox = new StsGroupBox("Geographically Limit the Wells Loaded");
    StsBooleanFieldBean projectBoundsBean = new StsBooleanFieldBean();
    StsFloatFieldBean haloBean = new StsFloatFieldBean();

    StsListFieldBean wellList = new StsListFieldBean();
    JScrollPane wellScrollPane = new JScrollPane();

    Object[] selectedWells = null;
    String currentWell = null;
    String[] wells = null;
    float halo = 0.0f;

    public StsWellHeadDefinitionPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsWellWizard)wizard;
        this.wizardStep = (StsWellHeadDefinition)wizardStep;

        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        xBean.initialize(wizard, "x", true, "X= ");
        yBean.initialize(wizard, "y", true, "Y= ");
        zBean.initialize(wizard, "z", true, "Datum:");
        haloBean.initialize(this, "halo", true, "Halo(ft/m):");
        kbApplyBean.initialize(wizard,"applyKb"," Apply Datum Correction");
        projectBoundsBean.initialize(this,"isProjectBound"," Limit load to Project Bounds.");
        wells = wizard.getSelectedWells();
        wellList.initialize(this,"selectedWell",null,wells);
        wellList.setSelectedIndex(0);
        currentWell = wells[0];
        wizard.rebuild();        
    }

    public void rebuildWellList()
    {
        wells = wizard.getSelectedWells();
        wellList.setListItems(wells);
        wellList.setSelectedIndex(0);
        currentWell = wells[0];
    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.BOTH;
        gbc.weighty = 1.0;
        gbc.anchor = gbc.NORTH;
        wellScrollPane.getViewport().add(wellList, null);
        add(wellScrollPane);

        typeBox.addEndRow(xBean);
        typeBox.addEndRow(yBean);
        typeBox.addEndRow(kbApplyBean);        
        typeBox.addEndRow(zBean);
        
        limitBox.addEndRow(projectBoundsBean);
        limitBox.addEndRow(haloBean);

        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0;
        add(typeBox);
        add(limitBox);
    }

    public void setSelectedWell(Object wells)
    {
        selectedWells = wellList.getSelectedObjects();
        currentWell = (String)wells;
        xBean.setValue(wizard.getX());
        yBean.setValue(wizard.getY());
        zBean.setValue(wizard.getZ());
        kbApplyBean.setValue(wizard.getApplyKb());
    }
    public Object getSelectedWell()
    {
        return currentWell;
    }

    public int[] getSelectedIndices()
    {
        if(selectedWells == null) return null;
        
        int[] indices = new int[selectedWells.length];
        for(int j=0; j<selectedWells.length; j++)
        {
            for(int i=0; i<wells.length; i++)
            {
                if((String)selectedWells[j] == wells[i])
                {
                    indices[j] = i;
                    break;
                }
            }
        }
        return indices;
    }
    public void setIsProjectBound(boolean val)
    {
        wizard.setIsProjectBound(val);
        rebuildWellList();
    }
    public boolean getIsProjectBound() { return wizard.getIsProjectBound();}
    public float getHalo() { return halo; }
    public void setHalo(float val)
    {
        halo = val;
        setIsProjectBound(getIsProjectBound());  // Recompute the wells in bounds if halo is changed.
    }
}