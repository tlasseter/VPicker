package com.Sts.Actions.Wizards.Collaboration;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsCollaborationSelectPanel extends JPanel implements ActionListener
{
	StsCollaborationEntry[] entries;
	StsCollaborationEntry currentEntry = null;
	private ButtonGroup zGroup = new ButtonGroup();
	private JTextArea sessionDescription = new JTextArea();
	private StsListFieldBean availableSessionsListBean = new StsListFieldBean();
	private JLabel sessionsLabel = new JLabel();
	private JLabel ipAddressLabel = new JLabel();
	public JCheckBox newSessionChk = new JCheckBox();
	private GridBagLayout gridBagLayout1 = new GridBagLayout();
	private JList participantsList = new JList();
//  JFormattedTextField ipAddress = new JFormattedTextField();
//  JLabel portLabel = new JLabel();
//  JFormattedTextField port = new JFormattedTextField();
  JLabel selectedSessionLabel = new JLabel();

	public StsCollaborationSelectPanel()
	{
		try
		{
			initialize();
			jbInit();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	void jbInit() throws Exception
	{
		sessionDescription.setFont(new java.awt.Font("Dialog", 0, 11));
		sessionDescription.setText("Session Description");
		this.setMinimumSize(new Dimension(0, 0));
		this.setPreferredSize(new Dimension(0, 0));
		this.setLayout(gridBagLayout1);
		sessionsLabel.setHorizontalAlignment(SwingConstants.LEFT);
		sessionsLabel.setText("Available Sessions");
		ipAddressLabel.setText("IP Address:");
		newSessionChk.setText("New Session");
//		ipAddress.setText("jFormattedTextField1");
//		portLabel.setText("Port:");
//    port.setText("jFormattedTextField1");
    selectedSessionLabel.setText("Selected session");
    this.add(sessionsLabel,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 0, 0));
		this.add(availableSessionsListBean,  new GridBagConstraints(0, 1, 1, 3, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 10, 10, 10), 0, 0));
		this.add(newSessionChk,  new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 0, 0));
		this.add(sessionDescription,   new GridBagConstraints(1, 1, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 10, 10, 10), 90, 82));
		this.add(participantsList,   new GridBagConstraints(1, 2, 2, 1, 1.0, 1.0
			,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 10, 10, 10), 90, 82));


/*
			this.add(ipAddressLabel,    new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 0, /));
		this.add(ipAddress,     new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 0, 0));
    this.add(portLabel,   new GridBagConstraints(1, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 0, 0));
    this.add(port,  new GridBagConstraints(2, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 0, 0));
*/
		this.add(selectedSessionLabel,  new GridBagConstraints(1, 0, 2, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 0, 0));
	}

	private void initialize()
	{
		StsCollaborationParticipant[] allParticipants = readAllParticipants();
		entries = readCollaborationEntries(allParticipants);
		availableSessionsListBean.initialize(this, "sessionSelected", null, entries);
	}

	static public StsCollaborationEntry[] readCollaborationEntries(StsCollaborationParticipant[] allParticipants)
	{
		String directory = System.getProperty("user.home") + File.separator + "S2SCache" + File.separator;
		String filename = "s2s.user.sessions";
		return StsCollaborationEntry.readFile(directory, filename, allParticipants);
	}

	static public StsCollaborationParticipant[] readAllParticipants()
	{
		String directory = System.getProperty("user.home") + File.separator + "S2SCache" + File.separator;
		String filename = "s2s.user.peers";
		return StsCollaborationParticipant.readFile(directory, filename);
	}

	public void setSessionSelected(Object entrySelected)
	{
		currentEntry = (StsCollaborationEntry)entrySelected;
		if(currentEntry == null) return;
		sessionDescription.setText(currentEntry.description);
		participantsList.setListData(currentEntry.participants);
		String username = System.getProperty("user.name");
		StsCollaborationParticipant participant = StsCollaborationParticipant.getParticipantWithUsername(currentEntry.participants, username);
		participantsList.setSelectedValue(participant, true);
		StsException.systemError("StsCollaborationSelectPanel.setSessionSelected() failed to match username.");
	}

	public Object getSessionSelected() { return currentEntry; }

	public void actionPerformed(ActionEvent e)
	{
	}
}
