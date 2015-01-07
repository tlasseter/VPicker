package com.Sts.Actions.Wizards.AnalogueCube;

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
	private StsSubVolume currentSubVolume;

	private StsSeismicVolume[] availableVolumes = null;
	private StsSeismicVolume volume;
	private String optionalName = "";
	private String cubeName = "analogueCube";
	private StsModel model = null;

	JLabel message = new JLabel("Attribute will need to be calculated");

	JList subVolumeList = new JList();
	DefaultListModel subVolumeListModel = new DefaultListModel();

	StsJPanel beanPanel = StsJPanel.addInsets();
	StsComboBoxFieldBean volumeListBean = new StsComboBoxFieldBean(this, "volume", "Target PostStack3d: ");
	StsStringFieldBean optionalNameBean = new StsStringFieldBean(this, "optionalName", "", true, "Optional Name: ");
	StsStringFieldBean cubeNameBean = new StsStringFieldBean(this.getClass(), "cubeName", "", false, "Cube Name: ");
	StsComboBoxFieldBean subVolumeListBean = new StsComboBoxFieldBean(this, "subVolume", "SubVolume Constraint: ");

	JRadioButton byteDataRadio = new JRadioButton();
	JRadioButton floatDataRadio = new JRadioButton();
	ButtonGroup dataFormatGroup = new ButtonGroup();

	GridBagLayout gridBagLayout1 = new GridBagLayout();
	GridBagLayout gridBagLayout2 = new GridBagLayout();
	JPanel radioPanel = new JPanel();
	GridBagLayout gridBagLayout3 = new GridBagLayout();

	public StsSelectVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = wizardStep;

		try
		{
			jbInit();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	public void initialize()
	{
		boolean found = false;
		int numVols = 0;

		model = wizard.getModel();
		availableVolumes = (StsSeismicVolume[]) model.getCastObjectList(StsSeismicVolume.class);
		volume = availableVolumes[0];
		volumeListBean.setListItems(availableVolumes);
		setCubeName(volume.stemname);

		StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
		StsSubVolume[] subVolumes = subVolumeClass.getSubVolumes();
		int nSubVolumes = subVolumes.length;
		StsSubVolume[] subVolumeList = new StsSubVolume[nSubVolumes+1];
		subVolumeList[0] = new StsSubVolume(false, "Entire PostStack3d");
		for(int n = 0; n < nSubVolumes; n++)
			subVolumeList[n+1] = (StsSubVolume)subVolumes[n];
		subVolumeListBean.setListItems(subVolumeList);
     }

	void jbInit() throws Exception
	{
		this.setLayout(gridBagLayout1);
		radioPanel.setBorder(BorderFactory.createEtchedBorder());
		radioPanel.setLayout(gridBagLayout3);
		byteDataRadio.setText("8-bit data");
		floatDataRadio.setText("32-bit float data");
		byteDataRadio.setSelected(true);
		this.add(beanPanel,
				 new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.BOTH,
										new Insets(5, 5, 5, 5), 0, 0));
		beanPanel.setBorder(BorderFactory.createEtchedBorder());
//		beanPanel.initializeLayout();
		beanPanel.addEndRow(volumeListBean);
		beanPanel.addEndRow(subVolumeListBean);
		beanPanel.addEndRow(optionalNameBean);
		beanPanel.addEndRow(cubeNameBean);
		message.setFont(new java.awt.Font("Dialog", 1, 11));
		message.setForeground(Color.BLUE);

		radioPanel.add(byteDataRadio,
					   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.BOTH,
											  new Insets(5, 5, 5, 5), 0, 0));
		radioPanel.add(floatDataRadio,
					   new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.BOTH,
											  new Insets(5, 5, 5, 5), 0, 0));

		beanPanel.add(radioPanel,
					  new GridBagConstraints(0, 4, 2, 1, 1.0, 1.0, GridBagConstraints.NORTHWEST,
											 GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));

		beanPanel.add(message,
					  new GridBagConstraints(0, 5, 2, 1, 1.0, 1.0, GridBagConstraints.NORTH, GridBagConstraints.CENTER,
											 new Insets(5, 20, 5, 5), 0, 0));

		dataFormatGroup.add(byteDataRadio);
		dataFormatGroup.add(floatDataRadio);
//		optionalNameBean.setEnabled(true);
//		optionalNameBean.setBackground(SystemColor.menu);
	}

	public void actionPerformed(ActionEvent e)
	{
        Object source = e.getSource();
    }

	public void setSubVolume(StsSubVolume subVolume)
	{

        if (currentSubVolume != null)
        {
            currentSubVolume.setIsVisible(false);
			currentSubVolume.setIsApplied(false);
        }

		if (subVolume.getName().equals("Entire PostStack3d"))
		{
			currentSubVolume = null;
			return;
		}
		currentSubVolume = subVolume;

		currentSubVolume.setIsVisible(true);
    }

	public StsSubVolume getSubVolume()
	{
		return currentSubVolume;
	}

	public void setVolume(StsSeismicVolume topVolume)
	{
		volume = topVolume;
		setCubeName(volume.stemname + optionalName + ".analog");
		boolean isDataFloat = topVolume.getIsDataFloat();
		floatDataRadio.setSelected(isDataFloat);
		floatDataRadio.setEnabled(isDataFloat);
		byteDataRadio.setEnabled(isDataFloat);
	}

	public StsSeismicVolume getVolume()
	{
		return volume;
	}

	public void setOptionalName(String name)
	{
		optionalNameBean.setText(name);
		if(name.length() > 0)
			optionalName = new String("." + name);
		else
			optionalName = "";
		setCubeName(volume.stemname + optionalName + ".analog");
	}

	public String getOptionalName()
	{
		return optionalNameBean.getText();
	}

	public void setCubeName(String name)
	{
		cubeNameBean.setText(name);
	}

	public String getCubeName()
	{
		return cubeNameBean.getText();
	}

	public boolean isDataFloat()
	{
		return floatDataRadio.isSelected();
	}
}
