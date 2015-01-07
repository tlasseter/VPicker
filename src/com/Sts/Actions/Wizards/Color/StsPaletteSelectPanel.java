package com.Sts.Actions.Wizards.Color;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
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

public class StsPaletteSelectPanel extends JPanel implements WindowListener, ActionListener, ListSelectionListener, ComponentListener
{
    private JLabel jLabel1 = new JLabel();
    private JPanel selectedPalettePanel = new JPanel();
    private JPanel selectionPanel = new JPanel();
    private JPanel transferPanel = new JPanel();
    private JButton addBtn = new JButton();
    private JButton addAllBtn = new JButton();
    private JButton removeBtn = new JButton();
    private JButton removeAllBtn = new JButton();

    private DefaultListModel selectedPalettesListModel = new DefaultListModel();
    private DefaultListModel availablePalettesListModel = new DefaultListModel();
    private JTextField currentDirectoryField = new JTextField();
    private JButton directoryBrowseButton = new JButton();

    private JFileChooser chooseDirectory = null;

    private JScrollPane availablePalettesScrollPane = new JScrollPane();
    private JList availablePalettesList = new JList();
    private JScrollPane selectedPalletteScrollPane = new JScrollPane();
    private JList selectedPalettesList = new JList();

    private StsPaletteWizard wizard;
    private StsPaletteSelect wizardStep;

    private String currentDirectory = null;
    public String currentPalette = null;
    public StsColor[] currentColors = null;

    public Color[] newColors = new Color[255];

    private StsAbstractFileSet availableFileSet;
    private StsAbstractFile[] availableFiles;
    private StsAbstractFile[] selectedFiles;
    private String[] selectedPalettes;

    protected StsSpectrum spectrum;
    public StsColorscalePanel colorscalePanel = null;

    private static String group = "palette";
    private static String format = "txt";

    String[] asciiFiles = null;
    String[] fileEndings = null;

    private GridBagLayout gridBagLayout1 = new GridBagLayout();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();
    private GridBagLayout gridBagLayout3 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();

    public StsPaletteSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsPaletteWizard)wizard;
        this.wizardStep = (StsPaletteSelect)wizardStep;

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
        availablePalettesList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);

        colorscalePanel = new StsColorscalePanel(newColors);
        selectedPalettePanel.add(colorscalePanel);
    }

    private void setAvailableFiles()
    {
      // Clear the List
      availablePalettesListModel.removeAllElements();

      // Get palettes from disk
      availableFileSet = StsFileSet.constructor(currentDirectory, new paletteFilter(group, format));
      String[] palFiles = availableFileSet.getFilenames();

      if(palFiles != null) {
        for(int i=0; i<palFiles.length; i++) {
           StsKeywordIO.parseAsciiFilename(palFiles[i]);
           availablePalettesListModel.addElement(StsKeywordIO.name);
        }
      }

   }
//
//  Get a list of all selected files...
//
    public boolean setSelectedFiles()
    {
      selectedFiles = new StsFile[selectedPalettesListModel.getSize()];
      selectedPalettes = new String[selectedPalettesListModel.getSize()];
      for (int i=0; i<selectedPalettesListModel.getSize(); i++)
      {
           selectedPalettes[i] = selectedPalettesListModel.getElementAt(i).toString();
           selectedFiles[i] = availableFileSet.getFile("palette.txt." + selectedPalettes[i]);
      }
      if(selectedFiles.length > 0)
          selectedPalettesList.setSelectedIndex(0);
      return true;
    }

    public String[] getSelectedPaletteNames() { return selectedPalettes; }
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
        selectedPalettesList.setMinimumSize(new Dimension(75, 200));
        selectedPalettesList.setModel(selectedPalettesListModel);
        availablePalettesList.setModel(availablePalettesListModel);
        addBtn.setText("Add >");
        addBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addBtn.setMaximumSize(new Dimension(100, 20));
        addBtn.setMinimumSize(new Dimension(100, 20));
        addBtn.setPreferredSize(new Dimension(100, 20));
        addBtn.setMargin(new Insets(0, 0, 0, 0));
        addBtn.addActionListener(this);
        availablePalettesList.addListSelectionListener(this);
        selectedPalettesList.addListSelectionListener(this);
        this.addComponentListener(this);
        availablePalettesScrollPane.setPreferredSize(new Dimension(150, 110));
        selectedPalletteScrollPane.setPreferredSize(new Dimension(100, 150));
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
    selectionPanel.add(availablePalettesScrollPane,  new GridBagConstraints(0, 1, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(9, 3, 9, 0), -16, 61));
    selectionPanel.add(selectedPalletteScrollPane,  new GridBagConstraints(3, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(9, 6, 9, 7), 34, 21));
    selectionPanel.add(directoryBrowseButton,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 0, 0), 0, 0));
    selectedPalletteScrollPane.getViewport().add(selectedPalettesList, null);
    availablePalettesScrollPane.getViewport().add(availablePalettesList, null);
    selectionPanel.add(transferPanel,     new GridBagConstraints(2, 1, 1, 1, 0.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(9, 0, 9, 0), 55, 20));

    }

    public void valueChanged(ListSelectionEvent e)
    {
       Object source = e.getSource();
       String name = null;

       if(source == availablePalettesList)
       {
           // Read File
           currentColors = readPalette(availableFileSet.getFile(availablePalettesList.getSelectedIndex()));
           currentPalette = availablePalettesListModel.getElementAt(availablePalettesList.getSelectedIndex()).toString();
       }
       else if(source == selectedPalettesList)
       {
           currentColors = readPalette(selectedFiles[selectedPalettesList.getSelectedIndex()]);
           currentPalette = selectedPalettesListModel.getElementAt(selectedPalettesList.getSelectedIndex()).toString();
       }
       computeNewColors();
       colorscalePanel.setColors(newColors);
    }

    private void computeNewColors() {
        int nColors = 255;

        int nColorIntervals = currentColors.length - 1;

        int nextColorIndex = 0;
        float[] nextColorRGBA = currentColors[0].getRGBA();
        float[] scale = new float[3];
        int nn = 0;
        for(int n = 0; n < nColorIntervals; n++)
        {
            int nColorsInInterval = currentColors[n+1].idx - currentColors[n].idx;
            int prevColorIndex = nextColorIndex;
            nextColorIndex += nColorsInInterval;
            if(nextColorIndex > nColors)
                nextColorIndex = nColors;
            float[] prevColorRGBA = nextColorRGBA;
            nextColorRGBA = currentColors[n+1].getRGBA();

            for(nn = prevColorIndex; nn < nextColorIndex; nn++)
            {
                if(nColorsInInterval <= 1) {
                    newColors[nn] = new Color(nextColorRGBA[0], nextColorRGBA[1],
                                              nextColorRGBA[2], nextColorRGBA[3]);
                    nextColorIndex = currentColors[n+1].idx;
                    nextColorRGBA = currentColors[n+1].getRGBA();
                    continue;
                } else {
                    float temp = (float) (nn - prevColorIndex) / (float) nColorsInInterval;
                    float[] interpolatedColor = StsMath.interpolate( prevColorRGBA, nextColorRGBA, temp);
                    newColors[nn] = new Color(interpolatedColor[0], interpolatedColor[1],
                                              interpolatedColor[2], interpolatedColor[3]);
                }
            }
        }
        for( ; nn < nColors; nn++)
        {
            newColors[nn] = new Color(nextColorRGBA[0], nextColorRGBA[1],
                                      nextColorRGBA[2], nextColorRGBA[3]);
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

        selectedPalettesList.removeListSelectionListener(this);
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
                    if(!StsYesNoDialog.questionValue(this,"Must select the directory that\n contains the Palette Files.\n\n Continue?"))
                        break;
                }
            }
//
//            if(chooseDirectory == null) initializeChooseDirectory();
//            {
//                chooseDirectory = new JFileChooser(currentDirectory);
//                chooseDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
//            }
//            chooseDirectory.showOpenDialog(null);
//            File newDirectory = chooseDirectory.getSelectedFile();
//            if(newDirectory != null) setCurrentDirectory(newDirectory.getAbsolutePath());
        }
        // Add selected volumes
        else if(source == addBtn)
        {
            selectedIndices = availablePalettesList.getSelectedIndices();
            for(i=0; i< selectedIndices.length; i++)
            {
                if(selectedPalettesListModel.indexOf(availablePalettesList.getModel().getElementAt(selectedIndices[i])) >= 0)
                    continue;
                selectedPalettesListModel.addElement(availablePalettesList.getModel().getElementAt(selectedIndices[i]));
            }
        }
        // Remove selected volume(s)
        else if (source == removeBtn)
        {
            selectedIndices = selectedPalettesList.getSelectedIndices();
            for(i= selectedIndices.length - 1; i>=0; i--)
                selectedPalettesListModel.removeElementAt(selectedIndices[i]);
        }
        // Remove all volumes
        else if (source == removeAllBtn)
        {
            selectedPalettesListModel.removeAllElements();
        }
        // Add all volumes

        else if (source == addAllBtn)
        {
            selectedPalettesListModel.removeAllElements();
            for (i=0; i<availablePalettesList.getModel().getSize(); i++)
            {
                selectedPalettesListModel.addElement(availablePalettesList.getModel().getElementAt(i));
            }
        }
        selectedPalettesList.addListSelectionListener(this);
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


    static final class paletteFilter implements FilenameFilter
    {
       String format = null;
       String group = null;
       int length;

       public paletteFilter(String grp, String fmt)
       {
           this.format = fmt;
           this.group = grp;
           length = group.length() + format.length() + 2;  // string length of "prefix.format."
       }

       public boolean accept(File dir, String name)
       {
           StsKeywordIO.parseAsciiFilename(name);

           String grp = StsKeywordIO.group;
           String fmt = StsKeywordIO.format;

           return StsKeywordIO.format.equals(format) && StsKeywordIO.group.equals(group);
        }

        public String getFilenameEnding(String filename)
        {
           return filename.substring(length);
        }
    }

    public StsColor[] readPalette(StsAbstractFile file)
    {
        StsAsciiFile asciiFile = null;
            int r, g, b, a, nKeys = 0;
            int idx = -1;
            int arraySizeInc = 256;
            StsColor color = null;
            StsColor[] colors = null;

            asciiFile = new StsAsciiFile(file);
            if(!asciiFile.openReadWithErrorMessage()) return null;

            String filename = file.getFilename();

            String[] tokens;
            try
            {
              while( (tokens = asciiFile.getTokens()) != null)
              {
                int nTokens = tokens.length;
                if(nTokens < 3)
                {
                  String inputLine = asciiFile.getLine();
                  StsMessageFiles.errorMessage("Insufficient entries for line: " + inputLine);
                }
                else if(nTokens == 4) {
                  r = Integer.parseInt(tokens[0]);
                  g = Integer.parseInt(tokens[1]);
                  b = Integer.parseInt(tokens[2]);
                  a = Integer.parseInt(tokens[3]);
                  color = new StsColor(nKeys, r, g, b, a);
                } else {
                  idx = Integer.parseInt(tokens[0]);
                  r = Integer.parseInt(tokens[1]);
                  g = Integer.parseInt(tokens[2]);
                  b = Integer.parseInt(tokens[3]);
                  a = Integer.parseInt(tokens[4]);
                  color = new StsColor(idx, r, g, b, a);
                }
 /*               if(colors == null)
                {
                    colors = new StsColor[1];
                    colors[0] = color;
                }
                else
 */
                colors = (StsColor[])StsMath.arrayAddElement(colors, color, nKeys, arraySizeInc);

                if(colors == null)
                {
                  StsException.systemError("Cannot create palette for: " + filename +
                                 " due to StsMath.arrayAddElement error.");
                  return null;
                }
                nKeys++;
              }
            } catch (java.io.IOException e) {
              ;
            }
            asciiFile.close();
            if(nKeys > 255)
            	nKeys = 255;
            colors = (StsColor[])StsMath.trimArray(colors, nKeys);

            return colors;
        }

        public void windowClosing(WindowEvent e) {
            spectrum.delete();
        }
        public void windowDeactivated(WindowEvent e) { }
        public void windowOpening(WindowEvent e) { }
        public void windowActivated(WindowEvent e) { }
        public void windowDeiconified(WindowEvent e) { }
        public void windowIconified(WindowEvent e) { }
        public void windowClosed(WindowEvent e) { }
        public void windowOpened(WindowEvent e) { }

}
