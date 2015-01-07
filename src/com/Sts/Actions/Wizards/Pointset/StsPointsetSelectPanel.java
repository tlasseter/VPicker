package com.Sts.Actions.Wizards.Pointset;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsPointsetSelectPanel extends JPanel implements WindowListener, ActionListener, ListSelectionListener, ComponentListener
{
    private JLabel jLabel1 = new JLabel();
    private JPanel selectedPalettePanel = new JPanel();
    private JPanel selectionPanel = new JPanel();
    private JPanel transferPanel = new JPanel();
    private JButton addBtn = new JButton();
    private JButton addAllBtn = new JButton();
    private JButton removeBtn = new JButton();
    private JButton removeAllBtn = new JButton();

    private DefaultListModel selectedPointsetsListModel = new DefaultListModel();
    private DefaultListModel availablePointsetsListModel = new DefaultListModel();
    private JTextField currentDirectoryField = new JTextField();
    private JButton directoryBrowseButton = new JButton();

    private JFileChooser chooseDirectory = null;

    private JScrollPane availablePointsetsScrollPane = new JScrollPane();
    private JList availablePointsetsList = new JList();
    private JScrollPane selectedPointsetScrollPane = new JScrollPane();
    private JList selectedPointsetsList = new JList();

    private StsPointsetWizard wizard;
    private StsPointsetSelect wizardStep;

    private String currentDirectory = null;
    public String currentPointset = null;

    private StsAbstractFileSet availableFileSet;
    private StsAbstractFile[] availableFiles;
    private StsAbstractFile[] selectedFiles;
    private String[] selectedPointsets;

    protected StsPointList pointset;

    String[] asciiFiles = null;
    String[] fileEndings = null;

    private GridBagLayout gridBagLayout1 = new GridBagLayout();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();
    private GridBagLayout gridBagLayout3 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();

    public StsPointsetSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsPointsetWizard)wizard;
        this.wizardStep = (StsPointsetSelect)wizardStep;

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
        StsModel model = wizard.getModel();
        if(model != null)
            currentDirectory = wizard.getModel().getProject().getRootDirString();
        else
            currentDirectory = System.getProperty("user.dirNo"); // standalone testing

        currentDirectoryField.setText(currentDirectory);

        setAvailableFiles();
        availablePointsetsList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    }

    private void setAvailableFiles()
    {
      // Clear the List
      availablePointsetsListModel.removeAllElements();

      // Get palettes from disk
      availableFileSet = StsFileSet.constructor(currentDirectory, new pointsetFilter("csv"));
      String[] palFiles = availableFileSet.getFilenames();

      if(palFiles != null)
      {
        for(int i=0; i<palFiles.length; i++)
        {
           String name = palFiles[i].substring(0,palFiles[i].lastIndexOf(".csv"));
           availablePointsetsListModel.addElement(name);
        }
      }

   }
//
//  Get a list of all selected files...
//
    public boolean setSelectedFiles()
    {
      selectedFiles = new StsFile[selectedPointsetsListModel.getSize()];
      selectedPointsets = new String[selectedPointsetsListModel.getSize()];
      for (int i=0; i<selectedPointsetsListModel.getSize(); i++)
      {
           selectedPointsets[i] = selectedPointsetsListModel.getElementAt(i).toString();
           selectedFiles[i] = availableFileSet.getFile(selectedPointsets[i] + ".csv");
      }
      if(selectedFiles.length > 0)
          selectedPointsetsList.setSelectedIndex(0);
      return true;
    }

    public String[] getSelectedPointsetNames() { return selectedPointsets; }
    public StsAbstractFile[] getSelectedFiles() { return selectedFiles; }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout4);
        this.setMinimumSize(new Dimension(600, 600));
        selectedPalettePanel.setLayout(gridBagLayout3);
        selectionPanel.setBorder(BorderFactory.createEtchedBorder());
        selectionPanel.setLayout(gridBagLayout2);
        selectedPalettePanel.setMaximumSize(new Dimension(30, 449));
        selectedPalettePanel.setMinimumSize(new Dimension(30, 449));
        selectedPalettePanel.setPreferredSize(new Dimension(30, 300));
        this.setMinimumSize(new Dimension(0, 0));
        this.setPreferredSize(new Dimension(0, 0));
        removeAllBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeAllBtn.setMaximumSize(new Dimension(100, 20));
        removeAllBtn.setMinimumSize(new Dimension(100, 20));
        removeAllBtn.setPreferredSize(new Dimension(100, 20));
        removeAllBtn.setMargin(new Insets(0, 0, 0, 0));
        removeAllBtn.setText("<< Remove All");
        removeAllBtn.addActionListener(this);
        removeBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeBtn.setMaximumSize(new Dimension(100, 20));
        removeBtn.setMinimumSize(new Dimension(100, 20));
        removeBtn.setPreferredSize(new Dimension(100, 20));
        removeBtn.setMargin(new Insets(0, 0, 0, 0));
        removeBtn.setText("< Remove");
        removeBtn.addActionListener(this);
        addAllBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addAllBtn.setMaximumSize(new Dimension(100, 20));
        addAllBtn.setMinimumSize(new Dimension(100, 20));
        addAllBtn.setPreferredSize(new Dimension(100, 20));
        addAllBtn.setMargin(new Insets(0, 0, 0, 0));
        addAllBtn.setText("Add All >>");
        addAllBtn.addActionListener(this);
        transferPanel.setLayout(gridBagLayout1);
        transferPanel.setMaximumSize(new Dimension(100, 200));
        transferPanel.setMinimumSize(new Dimension(25, 140));
        transferPanel.setPreferredSize(new Dimension(50, 150));
        currentDirectoryField.setPreferredSize(new Dimension(100, 20));
//        currentDirectoryField.setText("jTextField1");
        ImageIcon icon = StsIcon.createIcon("dir16x32.gif");
        directoryBrowseButton.setIcon(icon);
        directoryBrowseButton.setMargin(new Insets(0, 0, 0, 0));
        directoryBrowseButton.addActionListener(this);
        selectedPointsetsList.setMinimumSize(new Dimension(75, 200));
        selectedPointsetsList.setModel(selectedPointsetsListModel);
        availablePointsetsList.setModel(availablePointsetsListModel);
        addBtn.setText("Add >");
        addBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addBtn.setMaximumSize(new Dimension(100, 20));
        addBtn.setMinimumSize(new Dimension(100, 20));
        addBtn.setPreferredSize(new Dimension(100, 20));
        addBtn.setMargin(new Insets(0, 0, 0, 0));
        addBtn.addActionListener(this);
        availablePointsetsList.addListSelectionListener(this);
        selectedPointsetsList.addListSelectionListener(this);
        this.addComponentListener(this);
        availablePointsetsScrollPane.setPreferredSize(new Dimension(150, 110));
        selectedPointsetScrollPane.setPreferredSize(new Dimension(100, 150));
    this.add(selectionPanel,   new GridBagConstraints(0, 0, 6, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 4, 0), 0, 0));
    transferPanel.add(removeBtn,        new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 5));
    transferPanel.add(addAllBtn,      new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 5));
    transferPanel.add(addBtn,      new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 0, 1), 0, 5));
    transferPanel.add(removeAllBtn,         new GridBagConstraints(0, 3, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 5));
    this.add(selectedPalettePanel,     new GridBagConstraints(6, 0, 1, 1, 0.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(5, 0, 4, 6), 0, 0));
    selectionPanel.add(currentDirectoryField,  new GridBagConstraints(1, 0, 3, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(3, 13, 0, 7), 0, 3));
    selectionPanel.add(availablePointsetsScrollPane,  new GridBagConstraints(0, 1, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(9, 3, 9, 0), -16, 61));
    selectionPanel.add(selectedPointsetScrollPane,  new GridBagConstraints(3, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(9, 6, 9, 7), 34, 21));
    selectionPanel.add(directoryBrowseButton,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 0, 0), 0, 0));
    selectedPointsetScrollPane.getViewport().add(selectedPointsetsList, null);
    availablePointsetsScrollPane.getViewport().add(availablePointsetsList, null);
    selectionPanel.add(transferPanel,     new GridBagConstraints(2, 1, 1, 1, 0.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(9, 0, 9, 0), 55, 20));

    }

    public void valueChanged(ListSelectionEvent e)
    {
       Object source = e.getSource();
       String name = null;

       if(source == availablePointsetsList)
       {
           currentPointset = availablePointsetsListModel.getElementAt(availablePointsetsList.getSelectedIndex()).toString();
       }
       else if(source == selectedPointsetsList)
       {
           currentPointset = selectedPointsetsListModel.getElementAt(selectedPointsetsList.getSelectedIndex()).toString();
       }
    }

    public void componentHidden(ComponentEvent e) {}
    public void componentShown(ComponentEvent e) { }
    public void componentMoved(ComponentEvent e) {}
    public void componentResized(ComponentEvent e) { }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        int[] selectedIndices = null;

        Object source = e.getSource();

        selectedPointsetsList.removeListSelectionListener(this);
        // select new directory
        if(source == directoryBrowseButton)
        {
            if (chooseDirectory == null)
                initializeChooseDirectory();

            chooseDirectory = new JFileChooser(currentDirectory);
            chooseDirectory.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            chooseDirectory.setDialogTitle("Select or Enter Desired Directory and Press Open");
            chooseDirectory.setApproveButtonText("Open Directory");
            while(true)
            {
                chooseDirectory.showOpenDialog(null);
                File newDirectory = chooseDirectory.getSelectedFile();
                if (newDirectory.isDirectory())
                {
                    setCurrentDirectory(newDirectory.getAbsolutePath());
                    break;
                }
                else
                {
                    // File or nothing selected, strip off file name
                    String dirString = newDirectory.getPath();
                    newDirectory = new File(dirString.substring(0,dirString.lastIndexOf(File.separator)));
                    if(newDirectory.isDirectory())
                    {
                        setCurrentDirectory(newDirectory.getAbsolutePath());
                        break;
                    }
                    if(!StsYesNoDialog.questionValue(this,"Must select the directory that\n contains the Pointset Files.\n\n Continue?"))
                        break;
                }
            }
        }
        // Add selected volumes
        else if(source == addBtn)
        {
            selectedIndices = availablePointsetsList.getSelectedIndices();
            for(i=0; i< selectedIndices.length; i++)
            {
                if(selectedPointsetsListModel.indexOf(availablePointsetsList.getModel().getElementAt(selectedIndices[i])) >= 0)
                    continue;
                selectedPointsetsListModel.addElement(availablePointsetsList.getModel().getElementAt(selectedIndices[i]));
            }
        }
        // Remove selected volume(s)
        else if (source == removeBtn)
        {
            selectedIndices = selectedPointsetsList.getSelectedIndices();
            for(i= selectedIndices.length - 1; i>=0; i--)
                selectedPointsetsListModel.removeElementAt(selectedIndices[i]);
        }
        // Remove all volumes
        else if (source == removeAllBtn)
        {
            selectedPointsetsListModel.removeAllElements();
        }
        // Add all volumes

        else if (source == addAllBtn)
        {
            selectedPointsetsListModel.removeAllElements();
            for (i=0; i<availablePointsetsList.getModel().getSize(); i++)
            {
                selectedPointsetsListModel.addElement(availablePointsetsList.getModel().getElementAt(i));
            }
        }
        selectedPointsetsList.addListSelectionListener(this);
        // Create new selected file list
        setSelectedFiles();
    }

    private void initializeChooseDirectory()
    {
        chooseDirectory = new JFileChooser(currentDirectory);
        chooseDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    }

    private void setCurrentDirectory(String directory)
    {
        currentDirectory = directory;
        currentDirectoryField.setText(currentDirectory);
        setAvailableFiles();
    }

    public String getSelectedDirectory() { return currentDirectory; }


    static final class pointsetFilter implements FilenameFilter
    {
       String suffix = null;
       int length;

       public pointsetFilter(String sfx)
       {
           this.suffix = sfx;
           length = suffix.length() + 1;
       }

       public boolean accept(File dir, String fname)
       {
           StringTokenizer sTokens = new StringTokenizer(fname, ".");
           if(!sTokens.hasMoreTokens()) return false;
           String name = sTokens.nextToken();
           if(!sTokens.hasMoreTokens()) return false;
           String sfx = sTokens.nextToken();

           return sfx.equals(this.suffix);
        }

        public String getFilenameEnding(String filename)
        {
           return filename.substring(length);
        }
    }

    public void windowClosing(WindowEvent e) { }
    public void windowDeactivated(WindowEvent e) { }
    public void windowOpening(WindowEvent e) { }
    public void windowActivated(WindowEvent e) { }
    public void windowDeiconified(WindowEvent e) { }
    public void windowIconified(WindowEvent e) { }
    public void windowClosed(WindowEvent e) { }
    public void windowOpened(WindowEvent e) { }

}
