package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */
public class StsPoststackUnitsPanel extends StsGroupBox
{
	protected StsSeismicWizard wizard = null;
	protected StsComboBoxFieldBean domainComboBean;
	protected StsComboBoxFieldBean hUnitsBean;
	protected StsComboBoxFieldBean vUnitsBean;
//	private StsBooleanFieldBean nmoChk;

	// Z Domain
	public static final byte TIME = 1;
	public static final byte DEPTH = 2;

	public StsPoststackUnitsPanel(StsSeismicWizard wizard)
	{
		super("Units");
		this.wizard = wizard;
		constructBeans();
		constructPanel();
	}

	public void initialize()
	{
        setZDomainList();
 //       domainComboBean.getValueFromPanelObject();
//		hUnitsBean.getValueFromPanelObject();
//		vUnitsBean.getValueFromPanelObject();
//		nmoChk.setValueFromPanelObject();
	}

	protected void constructBeans()
	{
		domainComboBean = new StsComboBoxFieldBean(this, "zDomainString", "Domain:", StsParameters.TD_STRINGS);
        domainComboBean.setToolTipText("Select the domain of the input file (time or depth)");
		hUnitsBean = new StsComboBoxFieldBean(wizard, "horzUnits", "Horizontal:", StsParameters.DIST_STRINGS);
        hUnitsBean.setToolTipText("Select the horizontal units of the input file");
		vUnitsBean = new StsComboBoxFieldBean(wizard, "vertUnits", "Vertical:");
        vUnitsBean.setToolTipText("Select the vertical units of the input file");
	}

	protected void constructPanel()
	{
		gbc.fill = GridBagConstraints.BOTH;
		gbc.gridwidth = 3;

		domainComboBean.getComboBox().setEditable(false);
		hUnitsBean.getComboBox().setEditable(false);
		vUnitsBean.getComboBox().setEditable(false);

		addToRow(domainComboBean);
		addToRow(hUnitsBean);
		addToRow(vUnitsBean);
	}

	public void setZDomainString(String domain)
	{
		wizard.setZDomainString(domain);
		setZDomainList();
	}

	private void setZDomainList()
	{
        String zDomainString = wizard.getZDomainString();
        domainComboBean.setSelectedItem(zDomainString);
        if (wizard.getZDomainString() == StsParameters.TD_TIME_STRING)
			vUnitsBean.setListItems(StsParameters.TIME_STRINGS);
		else
			vUnitsBean.setListItems(StsParameters.DIST_STRINGS);
	}

	public String getZDomainString()
	{
		return wizard.getZDomainString();
	}
}
