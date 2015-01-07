package com.Sts.Actions.Wizards.VolumeFilter;

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

public class StsSelectVolumePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;
    private StsSubVolume[] availableSubVolumes = null;
    private StsSubVolume subVolume;

    private StsSeismicVolume[] availableVolumes = null;
    private StsSeismicVolume volume;
    private String cubeName = "Attribute";
    private StsModel model = null;

    StsJPanel beanPanel = StsJPanel.addInsets();
    StsComboBoxFieldBean volumeListBean = new StsComboBoxFieldBean(this, "volume", "Source Volume: ");
//    StsStringFieldBean cubeNameBean = new StsStringFieldBean(this.getClass(), "cubeName", "Attribute Name: ", true);
//    StsComboBoxFieldBean subVolumeListBean = new StsComboBoxFieldBean(this, "subVolumeString", "SubVolume Constraint: ", null);

//    JRadioButton byteDataRadio = new JRadioButton();
//    JRadioButton floatDataRadio = new JRadioButton();
//    ButtonGroup dataFormatGroup = new ButtonGroup();

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
 //   JPanel radioPanel = new JPanel();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsSelectVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
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
        availableVolumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
        volumeListBean.setListItems(availableVolumes);
        volumeListBean.setSelectedIndex(0);
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout1);
//        radioPanel.setBorder(BorderFactory.createEtchedBorder());
//        radioPanel.setLayout(gridBagLayout3);
//        byteDataRadio.setText("8-bit visualization format");
//        floatDataRadio.setText("32-bit float format");
//        byteDataRadio.setSelected(true);
        this.add(beanPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
        beanPanel.setBorder(BorderFactory.createEtchedBorder());
//        beanPanel.initializeLayout();
        beanPanel.gbc.fill = beanPanel.gbc.HORIZONTAL;
        beanPanel.add(volumeListBean);
//        beanPanel.addFieldBean(subVolumeListBean);
//        beanPanel.addFieldBean(cubeNameBean);

//        radioPanel.add(byteDataRadio, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
//            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
//        radioPanel.add(floatDataRadio, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
//            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
//        beanPanel.add(radioPanel, new GridBagConstraints(0, 3, 2, 1, 1.0, 1.0
//            ,GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));

//        dataFormatGroup.add(byteDataRadio);
//        dataFormatGroup.add(floatDataRadio);
    }

    public void actionPerformed(ActionEvent e) { }

    public StsSubVolume getSubVolume()
    {
        return subVolume;
    }

    public void setVolume(StsSeismicVolume topVolume)
    {
        volume = topVolume;
//		boolean isDataFloat = topVolume.getIsDataFloat();
//		floatDataRadio.setSelected(isDataFloat);
//		floatDataRadio.setEnabled(isDataFloat);
//		byteDataRadio.setEnabled(isDataFloat);
    }

    public StsSeismicVolume getVolume() { return volume; }

//    public void setCubeName(String name)
//    {
//        cubeNameBean.setText(name);
//    }
//    public String getCubeName() { return cubeNameBean.getText(); }

    public boolean isDataFloat()
    {
    	return volume.getIsDataFloat();
//        return floatDataRadio.isSelected();
    }
}
