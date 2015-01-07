package com.Sts.Actions.Wizards.Culture;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;

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

public class StsCultureSetSelectPanel extends JPanel implements ActionListener
{
    private JLabel jLabel1 = new JLabel();
    private JPanel selectedObjectInfoPanel = new JPanel();
    private JPanel selectionPanel = new JPanel();
    private JPanel transferPanel = new JPanel();
    private JButton addBtn = new JButton();
    private JButton addAllBtn = new JButton();
    private JButton removeBtn = new JButton();
    private JButton removeAllBtn = new JButton();

    private DefaultListModel selectedObjectsListModel = new DefaultListModel();
    private DefaultListModel availableObjectsListModel = new DefaultListModel();
    private JTextField currentDirectoryField = new JTextField();
    private JButton directoryBrowseButton = new JButton();

    private JFileChooser chooseDirectory = null;

    private JScrollPane availableObjectsScrollPane = new JScrollPane();
    private JList availableObjectsList = new JList();
    private JScrollPane selectedObjectsScrollPane = new JScrollPane();
    private JList selectedSurfacesList = new JList();

    private StsCultureWizard wizard;
    private StsCultureSetSelect wizardStep;

    private String currentDirectory = null;
    private StsAbstractFileSet availableFileSet;
    private StsAbstractFile[] availableFiles;
    private StsFile[] selectedFiles;
    private String[] selectedSurfaces = null;

    private GridBagLayout gridBagLayout1 = new GridBagLayout();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();
    JButton defineBtn = new JButton();

    JPanel jPanel1 = new JPanel();
    GridBagLayout gridBagLayout5 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();

    public StsCultureSetSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsCultureWizard)wizard;
        this.wizardStep = (StsCultureSetSelect)wizardStep;

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

        setAvailableFiles(currentDirectory);
        availableObjectsList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    }

    private void setAvailableFiles(String directory)
    {
        // Clear the List
        availableObjectsListModel.removeAllElements();
        if(!wizard.initializeCultureFilenameSet(directory)) return;
        Object[] objectNames = wizard.getObjectNames();
        for(int i=0; i<objectNames.length; i++)
            availableObjectsListModel.addElement(objectNames[i]);
    }

    public String[] getSelectedObjectNames()
    {
      String[] selectedObjectNames = new String[selectedObjectsListModel.getSize()];
      for (int i=0; i<selectedObjectsListModel.getSize(); i++)
      {
           selectedObjectNames[i] = selectedObjectsListModel.getElementAt(i).toString();
      }
      return selectedObjectNames;
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout4);
        this.setMinimumSize(new Dimension(600, 600));
        jLabel1.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel1.setText("Selected Culture Set Information");
        selectedObjectInfoPanel.setBorder(BorderFactory.createEtchedBorder());
        selectedObjectInfoPanel.setLayout(gridBagLayout3);
        selectionPanel.setBorder(BorderFactory.createEtchedBorder());
        selectionPanel.setLayout(gridBagLayout2);
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
        currentDirectoryField.setText("Directory");
        ImageIcon icon = StsIcon.createIcon("dir16x32.gif");
        directoryBrowseButton.setIcon(icon);
        directoryBrowseButton.setMargin(new Insets(0, 0, 0, 0));
        directoryBrowseButton.addActionListener(this);
        selectedSurfacesList.setMinimumSize(new Dimension(75, 200));
        selectedSurfacesList.setModel(selectedObjectsListModel);
        availableObjectsList.setModel(availableObjectsListModel);
        addBtn.setText("Add >");
        addBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addBtn.setMaximumSize(new Dimension(100, 20));
        addBtn.setMinimumSize(new Dimension(100, 20));
        addBtn.setPreferredSize(new Dimension(100, 20));
        addBtn.setMargin(new Insets(0, 0, 0, 0));
        addBtn.addActionListener(this);
        availableObjectsScrollPane.setPreferredSize(new Dimension(150, 110));
        selectedObjectsScrollPane.setPreferredSize(new Dimension(100, 150));
        defineBtn.setText("Define/View");
        defineBtn.addActionListener(this);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout5);

        this.add(selectionPanel,  new GridBagConstraints(0, 1, 5, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 5, 0, 5), -3, 40));
        transferPanel.add(removeBtn,          new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 4, 0, 1), 0, 5));
        transferPanel.add(addAllBtn,        new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 4, 0, 1), 0, 5));
        transferPanel.add(addBtn,       new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 0, 1), 0, 5));
        transferPanel.add(removeAllBtn,             new GridBagConstraints(0, 3, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 4, 0, 1), 0, 5));
        transferPanel.add(defineBtn,        new GridBagConstraints(0, 4, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 4, 0, 1), 10, 1));
        this.add(selectedObjectInfoPanel,  new GridBagConstraints(0, 0, 5, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 2, 5), 0, 2));
        selectionPanel.add(currentDirectoryField,  new GridBagConstraints(1, 0, 3, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(3, 13, 0, 7), 227, 3));
        selectionPanel.add(availableObjectsScrollPane,  new GridBagConstraints(0, 2, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(9, 3, 9, 0), -16, 61));
        selectionPanel.add(selectedObjectsScrollPane,  new GridBagConstraints(3, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(9, 6, 9, 7), 34, 21));
        selectionPanel.add(directoryBrowseButton,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 0, 0), 0, 0));
        selectionPanel.add(transferPanel, new GridBagConstraints(2, 2, 1, 1, 0.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(9, 0, 9, 0), 55, 20));
        selectedObjectsScrollPane.getViewport().add(selectedSurfacesList, null);

        availableObjectsScrollPane.getViewport().add(availableObjectsList, null);
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        int[] selectedIndices = null;

        Object source = e.getSource();

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
                if(newDirectory == null) return;
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
                    if(!StsYesNoDialog.questionValue(this,"Must select the directory that\n contains the Surface Files.\n\n Continue?"))
                        break;
                }
            }
        }
        // Add selected volumes
        else if(source == addBtn)
        {
            selectedIndices = availableObjectsList.getSelectedIndices();
            for(i=0; i< selectedIndices.length; i++)
            {
                if(selectedObjectsListModel.indexOf(availableObjectsList.getModel().getElementAt(selectedIndices[i])) == -1)
                    selectedObjectsListModel.addElement(availableObjectsList.getModel().getElementAt(selectedIndices[i]));
            }
        }
        // Remove selected volume(s)
        else if (source == removeBtn)
        {
            selectedIndices = selectedSurfacesList.getSelectedIndices();
            for(i= selectedIndices.length - 1; i>=0; i--)
                selectedObjectsListModel.removeElementAt(selectedIndices[i]);
        }
        // Remove all volumes
        else if (source == removeAllBtn)
        {
            selectedObjectsListModel.removeAllElements();
        }
        // Add all volumes
        else if (source == addAllBtn)
        {
            selectedObjectsListModel.removeAllElements();
            for (i=0; i<availableObjectsList.getModel().getSize(); i++)
            {
                selectedObjectsListModel.addElement(availableObjectsList.getModel().getElementAt(i));
            }
        }
        // Reset Available Files
//        else if (source == originalChk)
//            setAvailableFiles();

        else if(source == defineBtn)
        {
            if(getSelectedObjectNames().length == 0)
                new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "Need to select at least one culture object file from the available list\n" +
                                            " (left) and place it in the selected list (right).");
        }
    }

    private void checkEnableDefineBtn()
    {
        if(selectedObjectsListModel.getSize() > 0)
            defineBtn.setEnabled(true);
        else
            defineBtn.setEnabled(false);
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
        setAvailableFiles(directory);
    }

    public String getCurrentDirectory() { return currentDirectory; }

    final class filenameFilter implements FilenameFilter
    {
      String suffix = null;
      public filenameFilter(String suffix) { this.suffix = suffix; }
      public boolean accept(File dir, String name) { return name.endsWith(suffix); }
    }

    final class SurfaceFilenameFilter implements FilenameFilter
    {
        String group = null;
        String format = null;
        int length;

        public SurfaceFilenameFilter(String group, String format)
        {
              this.group = group;
              this.format = format;
              length = group.length() + format.length() + 2;  // string length of "prefix.format."
        }

        public boolean accept(File dir, String name)
        {
            StsKeywordIO.parseBinaryFilename(name);
            String keywordGroup = StsKeywordIO.group;
            String keywordFormat = StsKeywordIO.format;

            boolean groupOK = keywordGroup.equals(group);
            boolean formatOK = keywordFormat.equals(format);

            return StsKeywordIO.group.equals(group) && StsKeywordIO.format.equals(format);
        }

        public String getFilenameEnding(String filename)
        {
            return filename.substring(length);
        }
    }
}
