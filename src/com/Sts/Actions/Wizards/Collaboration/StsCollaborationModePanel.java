package com.Sts.Actions.Wizards.Collaboration;

import com.Sts.Collaboration.*;

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

public class StsCollaborationModePanel extends JPanel implements ActionListener
{
	private ButtonGroup zGroup = new ButtonGroup();
	private JRadioButton slaveModeRadio = new JRadioButton();
	private JRadioButton masterModeRadio = new JRadioButton();
	private JRadioButton multipleModeRadio = new JRadioButton();
	private ButtonGroup collabModeGroup = new ButtonGroup();
	private JTextArea collabModeDesc = new JTextArea();
	private GridBagLayout gridBagLayout1 = new GridBagLayout();

	public StsCollaborationModePanel()
	{
		try
		{
			jbInit();
			initialize();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	void jbInit() throws Exception
	{
		collabModeDesc.setBorder(BorderFactory.createEtchedBorder());
		collabModeDesc.setEditable(false);
		collabModeDesc.setText("Current user has the ability to control the viewoint and modify the " + "data model.");
		collabModeDesc.setLineWrap(true);
		collabModeDesc.setWrapStyleWord(true);
		masterModeRadio.setText("Leader");
		masterModeRadio.setVerticalAlignment(SwingConstants.BOTTOM);
		masterModeRadio.setVerticalTextPosition(SwingConstants.BOTTOM);
		masterModeRadio.setSelected(true);
		masterModeRadio.addActionListener(this);

		slaveModeRadio.setText("Follower");
		slaveModeRadio.addActionListener(this);

		this.setLayout(gridBagLayout1);
                /*
		multipleModeRadio.setText("Multiple");
		multipleModeRadio.setVerticalAlignment(SwingConstants.TOP);
		multipleModeRadio.setVerticalTextPosition(SwingConstants.TOP);
		multipleModeRadio.addActionListener(this);
*/
		collabModeGroup.add(masterModeRadio);
		collabModeGroup.add(slaveModeRadio);
//		collabModeGroup.add(multipleModeRadio);
		this.add(collabModeDesc,
				 new GridBagConstraints(1, 0, 1, 2, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(17, 12, 21, 22), 0, 70));
//		this.add(multipleModeRadio,
//				 new GridBagConstraints(0, 2, 1, 1, 0.0, 1.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(6, 27, 44, 0), 3, -2));
		this.add(masterModeRadio,
				 new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(35, 27, 0, 0), 5, 0));
		this.add(slaveModeRadio,
				 new GridBagConstraints(0, 1, 1, 1, 0.0, 1.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(7, 27, 0, 0), 15, -2));
	}

	private void initialize()
	{
		StsCollaboration collaboration = StsCollaboration.getCollaboration();
		if(collaboration == null) return;
		if(collaboration.isLeader())
			masterModeRadio.setSelected(true);
		else
			slaveModeRadio.setSelected(true);
	}

	public void actionPerformed(ActionEvent e)
	{
		StsCollaboration collaboration = StsCollaboration.getCollaboration();
		if(collaboration == null) return;

		if (e.getSource() == masterModeRadio)
		{
			collabModeDesc.setText("Current user has the ability to control the viewoint and modify the data model.");
			collabModeDesc.repaint();
			collaboration.promoted();
		}
		else if (e.getSource() == slaveModeRadio)
		{
			collabModeDesc.setText("Leader will control current users viewpoint.");
			collabModeDesc.repaint();
			collaboration.demoted();
		}
		else if (e.getSource() == multipleModeRadio)
		{
			collabModeDesc.setText("Current user and any other user in multiple or leader mode can have independent viewpoints");
			collabModeDesc.repaint();
		}
	}
}
