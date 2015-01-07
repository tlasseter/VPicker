package com.Sts.Actions.Wizards.Seismic;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

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

public class StsSEGYFormatSaveDialog extends JDialog
{
	private SEGYFormatPanel[] segyFormatPanels;

	private StsJPanel panel = StsJPanel.addInsets();
	private JPanel buttonPanel = new JPanel();
	private JButton cancelButton = null;
	private JButton okayButton = new JButton("OK");

	public StsSEGYFormatSaveDialog(Frame owner, StsSEGYFormat[] segyFormats)
	{
		super(owner, "Save SEGY Formats", true);
		try
		{
			setDefaultCloseOperation(DISPOSE_ON_CLOSE);
			jbInit(segyFormats);
			pack();
			StsToolkit.centerComponentOnScreen(this);
		}
		catch (Exception exception)
		{
			exception.printStackTrace();
		}
	}

	private void jbInit(StsSEGYFormat[] segyFormats) throws Exception
	{
		panel.addToRow(new JLabel("Modified   "));
		panel.addToRow(new JLabel("Save   "));
		panel.addToRow(new JLabel("SEGY Format Name   "));
		panel.addEndRow(new JLabel("Save As SEGY Format Name"));

		segyFormatPanels = new SEGYFormatPanel[segyFormats.length];
		for (int i = 0; i < segyFormats.length; i++)
		{
			segyFormatPanels[i] = new SEGYFormatPanel(segyFormats[i]);
		}
		buttonPanel.setLayout(new FlowLayout());
		panel.addEndRow(buttonPanel, 4, 1.0);
		buttonPanel.add(okayButton);
		cancelButton = new JButton("Cancel");
		buttonPanel.add(cancelButton);
		cancelButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				cancelButtonAction();
			}
		});
		okayButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				okayButtonAction();
			}
		});

		getContentPane().add(panel);
	}

	private void okayButtonAction()
	{
		for (int i = 0; i < segyFormatPanels.length; i++)
		{
			if (segyFormatPanels[i].save)
			{
				segyFormatPanels[i].segyFormat.serializeSilently(segyFormatPanels[i].getNewName());
			}
		}
		setVisible(false);
	}

	private void cancelButtonAction()
	{
		setVisible(false);
	}

	private class SEGYFormatPanel
	{
		public String name;
		public boolean save;
		public StsSEGYFormat segyFormat;

		JCheckBox saveCB = new JCheckBox();
		JTextField newName = new JTextField(); ;

		public SEGYFormatPanel(StsSEGYFormat segyFormat)
		{
			this.segyFormat = segyFormat;
			panel.gbc.fill = GridBagConstraints.HORIZONTAL;
			name = segyFormat.getName();
			save = segyFormat.isFormatChanged();
			JCheckBox changedCB = new JCheckBox();
			changedCB.setSelected(save);
			changedCB.setEnabled(false);
			panel.addToRow(changedCB);
			saveCB.setSelected(save);
			saveCB.setEnabled(true);
			panel.addToRow(saveCB);
			JLabel existingName = new JLabel(name);
			panel.addToRow(existingName);

			newName.setEditable(save);
			if (name.equalsIgnoreCase("SEG-Y"))
			{
				name = StsObject.getCurrentProject().getName();
			}
			newName.setText(name);
			panel.addEndRow(newName);

			saveCB.addActionListener(new java.awt.event.ActionListener()
			{
				public void actionPerformed(ActionEvent e)
				{
					setEnabled();
				}
			});
		}

		public String getNewName()
		{
			return newName.getText();
		}

		private void setEnabled()
		{
			boolean save = saveCB.isSelected();
			newName.setEditable(save);
		}
	}
}
