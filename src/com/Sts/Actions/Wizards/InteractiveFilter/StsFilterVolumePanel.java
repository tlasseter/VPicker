package com.Sts.Actions.Wizards.InteractiveFilter;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

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

public class StsFilterVolumePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    StsFilterVirtualVolume virtualVolume = null;

    private StsModel model = null;

    StsGroupBox[] filterPanels = null;
    StsComboBoxFieldBean convKernelCombo = new StsComboBoxFieldBean();
    StsComboBoxFieldBean rankSubTypeCombo = new StsComboBoxFieldBean();

    StsDoubleFieldBean radiusBean = new StsDoubleFieldBean();
    double filterRadius = 1;

    JPanel jPanel1 = new JPanel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    JPanel jPanel2 = new JPanel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    StsGroupBox typeBox = new StsGroupBox("Type");
    ButtonGroup typeGrp = new ButtonGroup();
    JRadioButton[] typeBtns = null;

    public StsFilterVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
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
    /*
        typeBtns = new JRadioButton[StsFilterVirtualVolume.FILTERS.length];
        filterPanels = new StsGroupBox[StsFilterVirtualVolume.FILTERS.length];
        for(int i=0; i<StsFilterVirtualVolume.FILTERS.length; i++)
        {
            typeBtns[i] = new JRadioButton(StsFilterVirtualVolume.FILTERS[i]);
            typeBtns[i].addActionListener(this);

            filterPanels[i] = new StsGroupBox(StsFilterVirtualVolume.FILTERS[i] + " Parameters");
            typeBox.gbc.anchor = typeBox.gbc.WEST;
            typeGrp.add(typeBtns[i]);
            typeBox.add(typeBtns[i]);
        }

        convKernelCombo.initialize(this,"kernelString","Kernel:",StsConvolve.KERNELS);
        rankSubTypeCombo.initialize(this,"subTypeString","Sub-Type:",StsRankFilters.RANKFILTERS);
        radiusBean.initialize(this,"filterRadius", 1, 20, "Radius:", true);

        this.add(typeBox,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 2, 4, 5), 0, 0));

        filterPanels[StsFilterVirtualVolume.CONVOLUTION].add(convKernelCombo);
        filterPanels[StsFilterVirtualVolume.RANK].add(rankSubTypeCombo);
        filterPanels[StsFilterVirtualVolume.RANK].add(radiusBean);

        typeBtns[StsFilterVirtualVolume.CONVOLUTION].setSelected(true);

        configureFilterPanels();
    */
    }

    private void configureFilterPanels()
    {
        /*
        for(int i=0; i<StsFilterVirtualVolume.FILTERS.length; i++)
            this.remove(filterPanels[i]);

        for(int i=0; i<StsFilterVirtualVolume.FILTERS.length; i++)
        {
            if(typeBtns[i].isSelected())
            {
                this.add(filterPanels[i],  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
                     ,GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(3, 2, 4, 5), 0, 0));
                validate();
                return;
            }
        }
        */
    }

    void jbInit() throws Exception
    {

    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source instanceof JRadioButton)
        {
            configureFilterPanels();
            virtualVolume = null;
        }
    }

    public void setKernelString(String kernel)
    {
        createEditDefinition();
    }
//    public String getKernelString() { StsConvolve.KERNELS[getKernel()]; }
    public byte getKernel()
    {
        return (byte)convKernelCombo.getSelectedIndex();
    }
    public void setSubTypeString(String subType)
    {
        createEditDefinition();
    }
 //   public String getSubTypeString() { return StsRankFilters.RANKFILTERS[getSubType()]; }
    public byte getSubType()
    {
        return (byte)rankSubTypeCombo.getSelectedIndex();
    }
    public int getFilterType()
    {
    /*
        for(int i=0; i<StsFilterVirtualVolume.FILTERS.length; i++)
        {
            if(typeBtns[i].isSelected())
                return i;
        }
    */
        return 0;
    }

    public double getFilterRadius() { return filterRadius; }
    public void setFilterRadius(double radius)
    {
        filterRadius = radius;
        createEditDefinition();
    }

    public boolean createEditDefinition()
    {
        /*
        if(virtualVolume == null)
        {
            StsSeismicBoundingBox[] volumes = new StsSeismicBoundingBox[] { StsPreStackLineSetClass.currentPreStackObject.velocityModel };
            if(getFilterType() == StsFilterVirtualVolume.CONVOLUTION)
                virtualVolume = new StsFilterVirtualVolume(volumes, "filterVolume", StsFilterVirtualVolume.CONVOLUTION, getKernel());
            else
                virtualVolume = new StsFilterVirtualVolume(volumes, "filterVolume", StsFilterVirtualVolume.RANK, getSubType(), getFilterRadius());
        }
        else
        {
            if(getFilterType() == StsFilterVirtualVolume.CONVOLUTION)
            {
                virtualVolume.setOperatorKernel(getKernel());
            }
            else
            {
                virtualVolume.setOperatorSubType(getSubType());
                virtualVolume.setOperatorSize(getFilterRadius());
            }
        }
        */
        // Refresh the current views with new filter volume
        System.out.println("Update the views with new definition");
        return true;
    }
}
