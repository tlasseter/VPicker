package com.Sts.UI;

import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsDirectorySelectionDialog extends JDialog
{
	// panel containing components
	StsJPanel panel = StsJPanel.addInsets();

	// large format title
	JLabel fileViewLbl = null;

	// resolution selection; if scaled float selected, sampleScale is enabled; text explanation displayed for each selection
	StsGroupBox directoryGroupBox = new StsGroupBox("Select output directory");
    StsButton directoryBrowseButton = new StsButton("dir16x32", "Browse for output directory.", this, "directoryBrowse");
	StsStringFieldBean currentDirectoryBean = new StsStringFieldBean(this, "currentDirectory", "Directory: ");

	// execute/cancel buttons
	StsGroupBox buttonGroupBox = new StsGroupBox();
	StsButton processButton = new StsButton("Process", "Start the export process.", this, "process");
	StsButton cancelButton = new StsButton("Cancel", "Cancel this operation.", this, "cancel");

	Frame frame;
    String title = "Export";
    String name = "";
    boolean process = true;
    private JFileChooser chooseDirectory = null;
    private String currentDirectory = null;

    public StsDirectorySelectionDialog(Frame frame, String title, boolean modal)
    {
        super(frame, title, modal);
        this.frame = frame;
        this.title = title;
        try
        {
            jbInit();
            this.setSize(400, 180);
            validate();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsDirectorySelectionDialog(Frame frame, String title, String name, boolean modal)
    {
        this(frame, title, modal);
        setName(name);
    }

    public StsDirectorySelectionDialog(Frame frame, String title, String directory, String filename, boolean modal)
    {
        this(frame, title, modal);
        setName(filename);
        setCurrentDirectory(directory);
    }

    private void jbInit() throws Exception
    {
        fileViewLbl = new JLabel(title, SwingConstants.CENTER);
		fileViewLbl.setFont(new java.awt.Font("Serif", 1, 14));

        panel.setBorder(BorderFactory.createEtchedBorder());
        this.getContentPane().add(panel);

        directoryGroupBox.gbc.weightx = 0.0;
        directoryGroupBox.gbc.fill = directoryGroupBox.gbc.HORIZONTAL;
        directoryGroupBox.addToRow(directoryBrowseButton);
        directoryGroupBox.gbc.weightx = 1.0;
        directoryGroupBox.addEndRow(currentDirectoryBean);

		buttonGroupBox.addToRow(processButton);
		buttonGroupBox.addEndRow(cancelButton);

        panel.gbc.fill = GridBagConstraints.HORIZONTAL;
	    panel.add(fileViewLbl);
		panel.add(directoryGroupBox);

	    buttonGroupBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		panel.add(buttonGroupBox);
    }

    private void initializeChooseDirectory()
    {
        chooseDirectory = new JFileChooser(currentDirectory);
        chooseDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    }

    public void setCurrentDirectory(String directory)
    {
        if(directory == null) return;
        if(!directory.endsWith(File.separator) && !directory.endsWith("\\"))
            directory = directory + File.separator;
        currentDirectory = directory;
        if(currentDirectoryBean != null)
            currentDirectoryBean.setValue(currentDirectory + name);
    }

    public String getCurrentDirectory()
    {
        return currentDirectory;
	}

    public void directoryBrowse()
    {
        int[] selectedIndices = null;

        if(chooseDirectory == null) initializeChooseDirectory();

        chooseDirectory = new JFileChooser(currentDirectory);
        chooseDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooseDirectory.setDialogTitle("Select or Enter Desired Directory and Press Accept");
        chooseDirectory.setApproveButtonText("Accept");
        while(true)
        {
            chooseDirectory.showOpenDialog(null);
            File newDirectory = chooseDirectory.getSelectedFile();
            if(newDirectory.isDirectory())
            {
                setCurrentDirectory(newDirectory.getAbsolutePath());
                break;
            }
        }
	}

	public void process()
	{
		process = true;
        setVisible(false);
    }

    public void cancel()
    {
		process = false;
        currentDirectory = null;
        setVisible(false);
    }

    public boolean getDoProcess() { return process;  }

}
