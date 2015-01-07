package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

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

public class StsCrossplotVVolumePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;

    private StsSeismicVolume[] seismicVolumes = null;
    private StsVirtualVolume[] virtualVolumes = null;
    private StsCrossplot[] crossplots = null;
    private int nSeismicVolumes = 0;
    private int nVirtualVolumes = 0;
    private int nCrossplots = 0;
    private int selectedXplotSeismicIndex = 0;
    private int selectedXplotIndex = 0;

    JPanel jPanel2 = new JPanel();
    JLabel seismicTwoLabel1 = new JLabel();
    JLabel seismicOneLabel1 = new JLabel();
    JComboBox xplotCombo = new JComboBox();
    JComboBox xplotSeismicCombo = new JComboBox();
    JLabel jLabel4 = new JLabel();
    ButtonGroup volumeType = new ButtonGroup();
    ButtonGroup crossplotGroup = new ButtonGroup();
    JRadioButton inclusiveBtn = new JRadioButton();
    JRadioButton exclusiveBtn = new JRadioButton();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
  GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsCrossplotVVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
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
        seismicVolumes = ((StsVirtualVolumeWizard)wizard).getSeismicVolumes();
        nSeismicVolumes = seismicVolumes.length;
        for(int v = 0; v < nSeismicVolumes; v++)
        {
            xplotSeismicCombo.addItem(seismicVolumes[v].getName());
        }
        virtualVolumes = ((StsVirtualVolumeWizard)wizard).getVirtualVolumes();
        nVirtualVolumes = virtualVolumes.length;
        for(int v = 0; v < nVirtualVolumes; v++)
        {
            xplotSeismicCombo.addItem(virtualVolumes[v].getName());
        }
        if((nSeismicVolumes + nVirtualVolumes) > 0)
        {
            setSelectedXplotSeismicVolume(0);
        }
        //
        // Add all existing crossplots
        //
        crossplots = ((StsVirtualVolumeWizard)wizard).getCrossplots();
        nCrossplots = crossplots.length;
        for(int v = 0; v < nCrossplots; v++)
        {
            xplotCombo.addItem(crossplots[v].getName());
        }
        if(nCrossplots > 0)
            setSelectedXplot(0);
        inclusiveBtn.setSelected(true);
    }

    // Set the seismic combobox to the current index
    private void setSelectedXplotSeismicVolume(int volumeIndex)
    {
        xplotSeismicCombo.setSelectedIndex(volumeIndex);
    }
    // Set the crossplot combobox to the current index
    private void setSelectedXplot(int xplotIndex)
    {
        xplotCombo.setSelectedIndex(xplotIndex);
    }

    // Get the current seismic index
    public int getXplotSeismicVolumeIndex()
    {
        return selectedXplotSeismicIndex;
    }
    // Get the current crossplot index
    public int getXplotIndex()
    {
        return selectedXplotIndex;
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);

        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout1);
        seismicTwoLabel1.setText("Crossplot");
        seismicTwoLabel1.setFont(new java.awt.Font("Dialog", 3, 11));
        seismicOneLabel1.setText("PostStack3d");
        seismicOneLabel1.setFont(new java.awt.Font("Dialog", 3, 11));
        xplotCombo.addActionListener(this);
        xplotCombo.setToolTipText("Select optional second seismic volume or specify scalar");
        xplotSeismicCombo.addActionListener(this);
        xplotSeismicCombo.setEnabled(true);
    xplotSeismicCombo.setToolTipText("Select first seismic volume");
        jLabel4.setText("AND");
        jLabel4.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel4.setFont(new java.awt.Font("Dialog", 1, 11));
        inclusiveBtn.setToolTipText("Include all common samples");
        inclusiveBtn.setText("Inclusive");
        exclusiveBtn.setToolTipText("Exclude all common samples");
        exclusiveBtn.setText("Exclusive");
        jPanel2.add(xplotSeismicCombo,   new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 4, 0, 0), 124, 4));
        jPanel2.add(seismicOneLabel1,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 5, 0, 0), 108, 0));
        jPanel2.add(xplotCombo,    new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 1, 0, 8), 124, 4));
        jPanel2.add(seismicTwoLabel1,   new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 4, 0, 8), 96, 0));
    jPanel2.add(jLabel4, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(7, 4, 7, 4), 8, -4));
    jPanel2.add(inclusiveBtn, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 0, 7, 3), 0, 0));
    jPanel2.add(exclusiveBtn, new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 25, 8, 57), 16, 0));

        crossplotGroup.add(inclusiveBtn);
        crossplotGroup.add(exclusiveBtn);
        this.add(jPanel2,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                                                  ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 3, 6, 5), 0, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == xplotSeismicCombo)
            selectedXplotSeismicIndex = xplotSeismicCombo.getSelectedIndex();
        else if(source == xplotCombo)
            selectedXplotIndex = xplotCombo.getSelectedIndex();
    }

    public StsSeismicVolume getSelectedXplotSeismicVolume()
    {
        if(selectedXplotSeismicIndex > seismicVolumes.length-1)
            return virtualVolumes[selectedXplotSeismicIndex - seismicVolumes.length];
        else
            return seismicVolumes[selectedXplotSeismicIndex];
    }
    public StsCrossplot getSelectedXplot()
    {
        return crossplots[selectedXplotIndex];
    }

    public boolean isInclusive()
    {
        if(inclusiveBtn.isSelected())
            return true;
        else
            return false;
    }

    public StsCrossplot getSelectedCrossplot()
    {
        return getSelectedXplot();
    }

}
