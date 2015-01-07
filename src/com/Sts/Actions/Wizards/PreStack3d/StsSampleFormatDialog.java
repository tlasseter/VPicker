package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.UI.Beans.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsSampleFormatDialog extends JDialog
{
	public static final int USE_DEFAULT = 0;
	public static final int USE_FILE = 1;
	public static final int SKIP_FILE = 2;

	private StsJPanel topPanel = StsJPanel.addInsets();
	private JPanel panel1 = new JPanel();
	private BorderLayout borderLayout1 = new BorderLayout();
	private JPanel jPanel1 = new JPanel();
	private JPanel jPanel2 = new JPanel();
	private JButton jButton1 = new JButton();
	private JTextArea jTextArea1 = new JTextArea();
	private BorderLayout borderLayout2 = new BorderLayout();
	private StsJPanel jPanel3 = StsJPanel.addInsets();
	private StsJPanel jPanel4 = StsJPanel.addInsets();
	private JRadioButton jRadioButton1 = new JRadioButton();
	private JRadioButton jRadioButton2 = new JRadioButton();
	private JRadioButton jRadioButton3 = new JRadioButton();
	private ButtonGroup buttonGroup1 = new ButtonGroup();
	private JCheckBox jCheckBox1 = new JCheckBox();
	private String filename;
	private String sampleFormat;
	private String defaultSampleFormat;

	public StsSampleFormatDialog(Frame owner, String title, boolean modal, String filename, String sampleFormat, String defaultSampleFormat)
	{
		super(owner, title, modal);
		try
		{
			this.filename = filename;
			this.sampleFormat = sampleFormat;
			this.defaultSampleFormat = defaultSampleFormat;
			jbInit();
		}
		catch (Exception exception)
		{
			exception.printStackTrace();
		}
	}

	public StsSampleFormatDialog(String filename, String sampleFormat, String defaultSampleFormat)
	{
		this(null, "Sample Format Conflict Dialog", true, filename, sampleFormat, defaultSampleFormat);
	}

	public boolean getDontAskAgain()
	{
		return jCheckBox1.isSelected();
	}

	public int getResponse()
	{
		return buttonGroup1.getSelection().getMnemonic();
	}

	private void jbInit() throws Exception
	{
		this.setSize(400, 300);
		this.setResizable(false);
		this.setLocationRelativeTo(null);
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);

		topPanel.gbc.insets = new Insets(8, 8, 2, 8);
		topPanel.gbc.fill = GridBagConstraints.BOTH;
		panel1.setLayout(borderLayout1);
		topPanel.add(panel1);
		getContentPane().add(topPanel);
		jButton1.setText("OK");
		jButton1.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				okButton_actionPerformed(e);
			}
		});
		jTextArea1.setBackground(SystemColor.activeCaptionBorder);
		jTextArea1.setFont(jPanel3.getFont());
		jTextArea1.setText(" The sample format in the binary header (" + sampleFormat + ") for file \"" + filename + "\" is\n not the same " +
						   "as the wizard sample format (" + defaultSampleFormat + ") .\n\n What would you like to do?");
		jTextArea1.setLineWrap(true);
		jTextArea1.setWrapStyleWord(true);
		jPanel2.setLayout(borderLayout2);
		jRadioButton1.setText("Use the wizard sample format (" + defaultSampleFormat + ") for ALL files");
		jRadioButton2.setText("Use the sample format from this file (" + sampleFormat + ") for ALL files");
		jRadioButton3.setText("Don't add this file (remove from list)");
		jCheckBox1.setText("Always use this option (for the remainder of this wizard session)");
		panel1.add(jPanel1, java.awt.BorderLayout.SOUTH);
		jPanel1.add(jButton1);
		panel1.add(jPanel2, java.awt.BorderLayout.CENTER);
		jPanel2.add(jTextArea1, java.awt.BorderLayout.NORTH);
		jPanel2.add(jPanel3, java.awt.BorderLayout.CENTER);
		jPanel3.gbc.fill = GridBagConstraints.NONE;
		jPanel4.gbc.fill = GridBagConstraints.NONE;
		jPanel4.gbc.anchor = GridBagConstraints.WEST;
		jPanel4.addEndRow(jRadioButton1);
		jPanel4.addEndRow(jRadioButton2);
		jPanel4.addEndRow(jRadioButton3);
		jPanel3.addEndRow(jPanel4);
		jPanel3.gbc.fill = GridBagConstraints.NONE;
		jPanel3.gbc.anchor = GridBagConstraints.WEST;
		jPanel3.addEndRow(jCheckBox1);
		jRadioButton1.setMnemonic(USE_DEFAULT);
		jRadioButton2.setMnemonic(USE_FILE);
		jRadioButton3.setMnemonic(SKIP_FILE);
		buttonGroup1.add(jRadioButton1);
		buttonGroup1.add(jRadioButton2);
		buttonGroup1.add(jRadioButton3);

		ChangeListener changeListener = new ChangeListener()
		{
			public void stateChanged(ChangeEvent e)
			{
				jRadioButton1StateChanged(e);
			}
		};
		jRadioButton1.addChangeListener(changeListener);
		jRadioButton2.addChangeListener(changeListener);
		jRadioButton3.addChangeListener(changeListener);
		jRadioButton1StateChanged(null);
	}

	public void jRadioButton1StateChanged(ChangeEvent e)
	{
		boolean okToClose = jRadioButton1.isSelected() || jRadioButton2.isSelected() || jRadioButton3.isSelected();
		jButton1.setEnabled(okToClose);
		if (okToClose)
		{
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		}
		else
		{
			setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		}
	}

	void okButton_actionPerformed(ActionEvent e)
	{
		dispose();
	}
}
