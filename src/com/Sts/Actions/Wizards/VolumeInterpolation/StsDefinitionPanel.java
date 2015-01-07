package com.Sts.Actions.Wizards.VolumeInterpolation;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefinitionPanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;

    private StsSeismicVolume[] seismicVolumes = null;
    private StsVirtualVolume[] virtualVolumes = null;

    private int nSeismicVolumes = 0;
    private int nVirtualVolumes = 0;

    private int selectedVolumeIndex = 0;

    StsCheckbox userNullsChkBox = new StsCheckbox("Volume Nulls","Interpolate through all values that equal user supplied volume null value.");
    StsCheckbox s2sNullsChkBox = new StsCheckbox("Padded Traces & Clipped Values","Interpolate an irregular survey to a rectangular area.");
    
    JComboBox volumeCombo = new JComboBox();
    JLabel volumeLabel = new JLabel();
    
    JPanel jPanel3 = new JPanel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsDefinitionPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
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
        model = wizard.getModel();
        //
        // Add all existing seismic volumes
        //
        seismicVolumes = ((StsVolumeInterpolationWizard)wizard).getSeismicVolumes();
        nSeismicVolumes = seismicVolumes.length;
        for(int v = 0; v < nSeismicVolumes; v++)
            volumeCombo.addItem(seismicVolumes[v].getName());
        virtualVolumes = ((StsVolumeInterpolationWizard)wizard).getVirtualVolumes();
        nVirtualVolumes = virtualVolumes.length;
        for(int v = 0; v < nVirtualVolumes; v++)
        	volumeCombo.addItem(virtualVolumes[v].getName());
        if((nSeismicVolumes + nVirtualVolumes) > 0)
            setSelectedVolume(0);
    }

    // Set the seismic combobox to the current index
    private void setSelectedVolume(int volumeIndex)
    {
    	volumeCombo.setSelectedIndex(volumeIndex);
    }

    // Get the current seismic index
    public int getSelectedVolumeIndex()
    {
        return selectedVolumeIndex;
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);
        volumeLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        volumeLabel.setText("Volume:");
        volumeCombo.setToolTipText("Select seismic volume to interpolate");
        volumeCombo.addActionListener(this);
        
        jPanel3.setBorder(BorderFactory.createEtchedBorder());
        jPanel3.setLayout(gridBagLayout1);

        jPanel3.add(volumeLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 5, 2), 0, 0));
        jPanel3.add(volumeCombo,     new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        jPanel3.add(userNullsChkBox, new GridBagConstraints(0, 1, 2, 1, 0.0, 0.0
            , GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 2), 0, 0));
        jPanel3.add(s2sNullsChkBox,     new GridBagConstraints(0, 2, 2, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 0), 0, 0));

        add(jPanel3);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == volumeCombo)
            selectedVolumeIndex = volumeCombo.getSelectedIndex();
    }

    public StsSeismicBoundingBox getSelectedVolume()
    {
        int idx = selectedVolumeIndex;
        idx = idx -1;
        if(idx == -1)
            return null;
        if(idx > seismicVolumes.length-1)
            return virtualVolumes[idx - seismicVolumes.length];
        else
            return seismicVolumes[idx];
    }

    public boolean interpolateUserNulls()
    {
    	return userNullsChkBox.isSelected();
    }
    
    public boolean interpolateS2SNulls()
    {
    	return s2sNullsChkBox.isSelected();
    }
}
