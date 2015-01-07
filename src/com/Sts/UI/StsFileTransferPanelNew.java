package com.Sts.UI;

import com.Sts.IO.*;
import com.Sts.IO.StsFilenameFilterFace;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.io.*;


/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2005</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author TJLasseter
 * @version c62e
 */
public class StsFileTransferPanelNew extends StsGroupBox
{
	/** instance listening to changes on this panel */
	private StsAbstractFileSet availableFileSet;
	private StsAbstractFile[] availableFiles;
	private String currentDirectory = null;
	private StsFilenameFilterFace filenameFilter;
	private StsFileTransferObjectFaceNew listener;

	private StsJPanel directoryGroupBox = StsJPanel.addInsets();
    private StsJPanel controlGroupBox = new StsGroupBox();

	private StsJPanel transferPanel = StsJPanel.addInsets();

	private StsButton directoryBrowseButton = new StsButton("dir16x32", "Browse for SEGY volume directory.", this, "directoryBrowse");
	private StsStringFieldBean currentDirectoryBean = null;
	private StsBooleanFieldBean reloadCheckbox;
    private StsBooleanFieldBean archiveItCheckbox;
    private StsBooleanFieldBean overrideCheckbox;

	private StsButton addBtn = new StsButton("Add >", "Add selected file to right side.", this, "addFiles");
	private StsButton removeBtn = new StsButton("< Remove", "Removed selected file from right side.", this, "removeFiles");
	private StsButton addAllBtn = new StsButton("Add all >", "Add all files to right side.", this, "addAllFiles");
	private StsButton removeAllBtn = new StsButton("< Remove all", "Remove all files from right side.", this, "removeAllFiles");

	private StsJPanel selectFilesPanel = StsJPanel.addInsets();
	private JFileChooser chooseDirectory = null;

	private DefaultListModel availableListModel = new DefaultListModel();
	private DefaultListModel selectedListModel = new DefaultListModel();
	private JList availableList = new JList(availableListModel);
	private JList selectedList = new JList(selectedListModel);

	private JScrollPane availableScrollPane = new JScrollPane();
	private JScrollPane selectedScrollPane = new JScrollPane();

    public StsFileTransferPanelNew(String currentDirectory, StsFilenameFilterFace filenameFilter, StsFileTransferObjectFaceNew listener,
                                int width, int height, boolean singleSelect)
    {
        this(currentDirectory, filenameFilter, listener, width, height);
        if(singleSelect)
        {
            addAllBtn.setEnabled(false);
            removeAllBtn.setEnabled(false);
            availableList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        }
	}

	public StsFileTransferPanelNew(String currentDirectory, StsFilenameFilterFace filenameFilter, StsFileTransferObjectFaceNew listener, int width, int height)
	{
		this(currentDirectory, filenameFilter, listener);
		int buttonPanelWidth = 100;
		transferPanel.setPreferredSize(new Dimension(buttonPanelWidth, height));
		int panelWidth = (width-buttonPanelWidth)/2;
		availableScrollPane.getViewport().setPreferredSize(new Dimension(panelWidth, height));
		selectedScrollPane.getViewport().setPreferredSize(new Dimension(panelWidth, height));
	}

	public StsFileTransferPanelNew(String currentDirectory, StsFilenameFilterFace filenameFilter, StsFileTransferObjectFaceNew listener)
	{
		super("File Transfer Selection");
        currentDirectoryBean = new StsStringFieldBean(this, "currentDirectory", "Directory: ");

		this.filenameFilter = filenameFilter;
		setCurrentDirectory(currentDirectory);
		this.listener = listener;

        if(listener.hasDirectorySelection())
        {
            directoryGroupBox.gbc.fill = GridBagConstraints.NONE;
            directoryGroupBox.gbc.anchor = GridBagConstraints.WEST;
            directoryGroupBox.gbc.weightx = 0.0;
            directoryGroupBox.gbc.weighty = 0.0;
            directoryGroupBox.addToRow(directoryBrowseButton);

            directoryGroupBox.gbc.weightx = 1.0;
            directoryGroupBox.gbc.fill = GridBagConstraints.HORIZONTAL;
            directoryGroupBox.gbc.anchor = GridBagConstraints.EAST;
            currentDirectoryBean.setColumns(30);
            directoryGroupBox.addEndRow(currentDirectoryBean);

            controlGroupBox.gbc.anchor = GridBagConstraints.WEST;
            boolean addControl = false;
            if (listener.hasReloadButton())
            {
                addControl = true;
                reloadCheckbox = new StsBooleanFieldBean(listener, "reload", listener.getReload(), "Reload Original Files");
                if(listener.hasArchiveItButton() || listener.hasOverrideButton())
                    controlGroupBox.addToRow(reloadCheckbox);
                else
                    controlGroupBox.addEndRow(reloadCheckbox);
            }
            if (listener.hasArchiveItButton())
            {
                addControl = true;
                archiveItCheckbox = new StsBooleanFieldBean(listener, "archiveIt", listener.getArchiveIt(), "Archive Files");
                if(listener.hasOverrideButton())
                    controlGroupBox.addToRow(archiveItCheckbox);
                else
                    controlGroupBox.addEndRow(archiveItCheckbox);
            }
            if (listener.hasOverrideButton())
            {
                addControl = true;
                overrideCheckbox = new StsBooleanFieldBean(listener, "overrideFilter", listener.getOverrideFilter(), "Override File Filter");
                if(listener.hasOverrideButton())
                    controlGroupBox.addToRow(overrideCheckbox);
                else
                    controlGroupBox.addEndRow(overrideCheckbox);
            }
            if(addControl)
            {
                directoryGroupBox.gbc.gridwidth = 3;
                directoryGroupBox.gbc.anchor = GridBagConstraints.WEST;
                directoryGroupBox.gbc.fill = GridBagConstraints.HORIZONTAL;
                directoryGroupBox.addEndRow(controlGroupBox);
            }
            gbc.anchor = GridBagConstraints.NORTH;
            gbc.fill = GridBagConstraints.HORIZONTAL;
            gbc.weighty = 0.0;
            addEndRow(directoryGroupBox);
        }
        transferPanel.gbc.weightx = 0.0;
		transferPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		transferPanel.add(addBtn);
		transferPanel.add(addAllBtn);
		transferPanel.add(removeBtn);
		transferPanel.add(removeAllBtn);

		availableScrollPane.getViewport().add(availableList, null);
		availableScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		availableScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

		selectedScrollPane.getViewport().add(selectedList, null);
		selectedScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		selectedScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

		selectFilesPanel.gbc.fill = GridBagConstraints.BOTH;
        selectFilesPanel.gbc.weighty = 1.0;
        selectFilesPanel.gbc.weightx = 0.5;
		selectFilesPanel.addToRow(availableScrollPane);
        selectFilesPanel.gbc.fill = GridBagConstraints.NONE;
        selectFilesPanel.gbc.weightx = 0.0;
		selectFilesPanel.addToRow(transferPanel);
        selectFilesPanel.gbc.fill = GridBagConstraints.BOTH;
        selectFilesPanel.gbc.weightx = 0.5;
		selectFilesPanel.addToRow(selectedScrollPane);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.weighty = 1.0;
		addEndRow(selectFilesPanel);

		availableList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		availableList.addListSelectionListener
		(
		    new ListSelectionListener()
	        {
				public void valueChanged(ListSelectionEvent e)
				{
					availableVolSelected(e);
				}
			}
		);
		selectedList.addListSelectionListener
		(
		    new ListSelectionListener()
	        {
				public void valueChanged(ListSelectionEvent e)
				{
					selectedVolSelected(e);
				}
			}
		);
	}
    public void removeFile(Object file)
	{
	    availableListModel.addElement(filenameFilter.getFilenameName(((StsAbstractFile)file).filename));
		selectedListModel.removeElement(filenameFilter.getFilenameName(((StsAbstractFile)file).filename));
    }
	public void availableVolSelected(ListSelectionEvent e)
	{
		Object source = e.getSource();
		if(!(source instanceof JList)) return;
        String selectedString = (String)availableList.getSelectedValue();
		if(selectedString == null) return;
		listener.availableFileSelected((StsFile)availableFileSet.getFileFromName(selectedString, filenameFilter));
	}

	public void selectedVolSelected(ListSelectionEvent e)
	{
		Object source = e.getSource();
		if(!(source instanceof JList)) return;
        String selectedString = (String)selectedList.getSelectedValue();
		if(selectedString == null) return;
		listener.fileSelected((StsFile)availableFileSet.getFileFromName(selectedString, filenameFilter));
	}

	public void selectSingleVolProgrammatically(String selectedString)
	{
		selectedList.setSelectedIndex(selectedListModel.indexOf(selectedString));
	}

	public void selectSingleVolProgrammatically(int index)
	{
		selectedList.setSelectedIndex(index);
	}

	public int getSelectedCount()
	{
		return selectedListModel.size();
	}

    public int[] getSelectedIndices() { return selectedList.getSelectedIndices(); }

    public void setFilenameFilter(StsFilenameFilterFace filter)
    {
        filenameFilter = filter;
    }

    public void setFilenameFilter(StsFilenameFilter filter)
    {
        filenameFilter = filter;
    }

	public void setAvailableFiles()
	{
        availableListModel.clear();
		availableFileSet = StsFileSet.constructor(currentDirectory, filenameFilter);
		availableFileSet.sort();
		availableFiles = availableFileSet.getFiles();
		for(int n = 0; n < availableFiles.length; n++)
        {
            availableListModel.addElement(filenameFilter.getFilenameName(availableFiles[n].getFilename()));
        }
    }

	public void directoryBrowse()
	{
		int[] selectedIndices = null;

		if(chooseDirectory == null) initializeChooseDirectory();

		chooseDirectory = new JFileChooser(currentDirectory);
		chooseDirectory.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		chooseDirectory.setDialogTitle("Select or Enter Desired Directory and Press Open");
		chooseDirectory.setApproveButtonText("Open Directory");
		while(true)
		{
			int retVal = chooseDirectory.showOpenDialog(null);
			if(retVal != chooseDirectory.APPROVE_OPTION)
				break;
			File newDirectory = chooseDirectory.getSelectedFile();
			if(newDirectory == null) return;
			if(newDirectory.isDirectory())
			{
				setCurrentDirectory(newDirectory.getAbsolutePath());
				break;
			}
			else
			{
				// File or nothing selected, strip off file name
				String dirString = newDirectory.getPath();
				newDirectory = new File(dirString.substring(0, dirString.lastIndexOf(File.separator)));
				if(newDirectory.isDirectory())
				{
					setCurrentDirectory(newDirectory.getAbsolutePath());
					break;
				}
				if(!StsYesNoDialog.questionValue(this,
											 "Must select the directory that\n contains the SegY Files.\n\n Continue?"))
					break;
			}
		}
	}

	private void initializeChooseDirectory()
	{
		chooseDirectory = new JFileChooser(currentDirectory);
		chooseDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	}

	public void setCurrentDirectory(String directory)
	{
        if(directory == null)
            return;
		currentDirectory = directory;
		if(currentDirectoryBean != null)
            currentDirectoryBean.setValue(currentDirectory);
		setAvailableFiles();
	}

	public String getCurrentDirectory()
	{
		return currentDirectory;
	}

	public void addFiles()
	{
		int[] selectedIndices = availableList.getSelectedIndices();
        if(selectedIndices.length == 0)
        {
            new StsMessage(this, StsMessage.INFO,  "No files selected.");
            return;
        }
        String[] addedStrings = new String[selectedIndices.length];
        StsFile[] addedFiles = new StsFile[selectedIndices.length];
		for (int i = 0; i < selectedIndices.length; i++)
			addedStrings[i] = (String)availableListModel.getElementAt(selectedIndices[i]);
		for (int i = 0; i < selectedIndices.length; i++)
		{
			selectedListModel.addElement(addedStrings[i]);
			availableListModel.removeElement(addedStrings[i]);
			addedFiles[i] = (StsFile)availableFileSet.getFileFromName(addedStrings[i], filenameFilter);
		}
		listener.addFiles(addedFiles);
	}

    public void moveToAvailableList()
    {
        int[] selectedIndices = selectedList.getSelectedIndices();
        String[] removedStrings = new String[selectedIndices.length];
        for (int i = 0; i < selectedIndices.length; i++)
            removedStrings[i] = (String)selectedListModel.getElementAt(selectedIndices[i]);
        for(int i = 0; i < selectedIndices.length; i++)
        {
            selectedListModel.removeElement(removedStrings[i]);
            availableListModel.addElement(removedStrings[i]);
        }
        if (selectedListModel.size() > 0) selectedList.setSelectedIndex(0);
    }

    public void removeFiles()
    {
        int[] selectedIndices = selectedList.getSelectedIndices();
        StsFile[] removedFiles = new StsFile[selectedIndices.length];
        String[] removedStrings = new String[selectedIndices.length];
        for (int i = 0; i < selectedIndices.length; i++)
            removedStrings[i] = (String)selectedListModel.getElementAt(selectedIndices[i]);
        for(int i = 0; i < selectedIndices.length; i++)
        {
            selectedListModel.removeElement(removedStrings[i]);
            availableListModel.addElement(removedStrings[i]);
            removedFiles[i] = (StsFile)availableFileSet.getFileFromName(removedStrings[i], filenameFilter);
        }
        if (selectedListModel.size() > 0) selectedList.setSelectedIndex(0);
        listener.removeFiles(removedFiles);
    }

	public void addAllFiles()
	{
		int nFiles = availableListModel.size();
        Object[] availableFiles = availableListModel.toArray();
        StsAbstractFile[] addedFiles  = (StsAbstractFile[])availableFileSet.getFiles();
        availableListModel.clear();
		listener.addFiles(addedFiles);
        for(int n = 0; n < nFiles; n++)
            selectedListModel.addElement(availableFiles[n]);
	}

	public void removeAllFiles()
	{
		int nFiles = selectedListModel.size();
		Object[] selectedFiles = selectedListModel.toArray();
		for(int i = 0; i < nFiles; i++)
		{
			String removedFile = (String)selectedFiles[i];
			availableListModel.addElement(removedFile);
		}
		selectedListModel.clear();
		listener.removeAllFiles();
	}

	public void removeSelectedFiles()
	{
		selectedListModel.clear();
		listener.removeAllFiles();
	}

	public void addAvailableFile(String fileName)
	{
		availableListModel.addElement(filenameFilter.getFilenameName(fileName));
	}

    public void refreshAvailableList()
    {
        setAvailableFiles();
    }

	public static void main(String[] args)
	{
		String[] filterStrings = new String[] {"sgy", "segy", "Segy", "SegY"};
		StsFilenameEndingFilter filter = new StsFilenameEndingFilter(filterStrings);
		String currentDirectory = System.getProperty("user.dirNo");
		TransferPanelTest panelTest = new TransferPanelTest();
		try{ UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel()); } catch(Exception e) { }
		StsFileTransferPanel panel = new StsFileTransferPanel(currentDirectory, filter, panelTest);
//		panel.setForeground(Color.gray);
//		panel.setFont(new java.awt.Font("Dialog", 1, 11));
		com.Sts.Utilities.StsToolkit.createDialog(panel);
	}
}