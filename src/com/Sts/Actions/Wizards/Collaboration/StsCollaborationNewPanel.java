package com.Sts.Actions.Wizards.Collaboration;

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

public class StsCollaborationNewPanel extends StsJPanel
{
	private StsStringFieldBean nameBean;
	private StsStringFieldBean addressBean;
	private StsIntFieldBean portBean;
	private StsStringFieldBean passwordBean;
	private JEditorPane descriptionText;

	public StsCollaborationNewPanel(StsCollaborationWizard wizard)
	{
		try
		{
			constructBeans(wizard);
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	void constructBeans(StsCollaborationWizard wizard) throws Exception
	{
		nameBean = new StsStringFieldBean(wizard, "name", "Name:");
		addressBean = new StsStringFieldBean(wizard, "ipAddress", "IP address:");
		portBean = new StsIntFieldBean(wizard, "port", true, "Port:");
		passwordBean = new StsStringFieldBean(wizard, "password", "Password:");
		JLabel descriptionLabel = new JLabel("Description:");
		descriptionLabel.setHorizontalAlignment(SwingConstants.RIGHT);
		JEditorPane descriptionText = new JEditorPane();
        gbc.fill = gbc.HORIZONTAL;
		add(nameBean);
		add(addressBean);
		add(portBean);
		add(passwordBean);
		this.gbc.gridwidth = 2;
		addToRow(descriptionLabel);
		addEndRow(descriptionText);
	}
/*
	public StsCollaborationEntry getNewEntry()
	{
		String name = nameText.getText();
		String ipAddress = ipAddressText.getText();
		int port = Integer.parseInt(portText.getText());
		String description = descriptionText.getText();
		[] participants = new String[] { "tom", "dick", "harry" };
		return new StsCollaborationEntry(name, description, participants);
	}
*/
}
