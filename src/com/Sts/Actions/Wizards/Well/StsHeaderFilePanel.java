package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

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

public class StsHeaderFilePanel extends StsJPanel implements ActionListener, ChangeListener {
    private StsWellWizard wizard;
    private StsHeaderFile wizardStep;
    
    private JFileChooser chooseDirectory = null;
    StsGroupBox fileSelectionBox = new StsGroupBox();
    JLabel jLabel4 = new JLabel("Selected File:");
    JButton dirBtn = new JButton();
    JTextField filenameTxt = new JTextField();

    JPanel jPanel1 = new JPanel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JScrollPane jScrollPane1 = new JScrollPane();
    JLabel fileViewLbl = new JLabel();
    JTextArea fileTextArea = new JTextArea();
    Font font = new java.awt.Font("Serif", 1, 10);
    int height = 500;
    int width = 500;
    private StsGroupBox colDefPanel = new StsGroupBox("# Header Rows");
    int numCols = 0;
    int numSkipped = 0;
    JLabel jLabel2 = new JLabel();
    JSpinner numSkippedSpin = new JSpinner();
    private SpinnerModel numSkippedSpinModel = null;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    private String currentFile = null;

    public StsHeaderFilePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsWellWizard)wizard;
        this.wizardStep = (StsHeaderFile)wizardStep;
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
        numSkippedSpin.setValue(new Integer(0));
        setValues();
    }

    void jbInit() throws Exception
    {
         this.setLayout(gridBagLayout1);
        jPanel1.setLayout(gridBagLayout2);

        numSkippedSpinModel = new SpinnerNumberModel(0, 0, 200, 1);
        
        dirBtn.setMargin(new Insets(2, 2, 2, 2));
        dirBtn.setText("Select...");
        dirBtn.addActionListener(this);
        filenameTxt.setText("");
        fileSelectionBox.addToRow(jLabel4,1,0.0);
        fileSelectionBox.gbc.fill = gbc.HORIZONTAL;
        fileSelectionBox.addToRow(filenameTxt,2,1.0);
        fileSelectionBox.gbc.fill = gbc.NONE;
        fileSelectionBox.addEndRow(dirBtn,1,0.0);

        fileViewLbl.setFont(new java.awt.Font("Serif", 1, 12));
        fileViewLbl.setHorizontalAlignment(SwingConstants.CENTER);
        fileViewLbl.setText("File View");
        fileTextArea.setBackground(Color.lightGray);
        fileTextArea.setEditable(false);
        fileTextArea.setText("First 100 lines of file will be isVisible here on selection...");
        fileViewLbl.setFont(font);
        fileTextArea.setLineWrap(false);
        fileTextArea.setFont(new Font("Monospaced", 1, 10));
        fileTextArea.setBorder(BorderFactory.createEtchedBorder());
        fileTextArea.setDoubleBuffered(true);
        fileTextArea.setCaretColor(Color.red);
        fileTextArea.setCaretPosition(0);
        jScrollPane1.setPreferredSize(new Dimension(width,height));

        jScrollPane1.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        jScrollPane1.getViewport().add(fileTextArea);
        jPanel1.add(jScrollPane1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(8, 6, 9, 7), 0, 0));

        numSkippedSpin.addChangeListener(this);
        numSkippedSpin.setModel(numSkippedSpinModel);

        colDefPanel.addEndRow(numSkippedSpin);

        gbc.fill = gbc.HORIZONTAL;
        gbc.gridheight = 1;
        gbc.weighty = 0.0;
        gbc.anchor = gbc.NORTH;
        this.addEndRow(fileSelectionBox,3,1.0);

        gbc.fill = gbc.BOTH;
        gbc.gridheight = 5;
        gbc.weighty = 1.0;
        this.addToRow(colDefPanel,1,0.0);
        this.addEndRow(jPanel1,2,1.0);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == dirBtn)
        {
            directoryBrowse();
            if(currentFile != null)
            {
                setValues();
            }
        }
    }

    public void directoryBrowse()
    {
        String currentDirectory = wizard.getModel().getProject().getRootDirString();
        chooseDirectory = new JFileChooser(currentDirectory);
        chooseDirectory.setFileSelectionMode(JFileChooser.FILES_ONLY);
        chooseDirectory.setDialogTitle("Select or Enter Desired File and Press Open");
        chooseDirectory.setApproveButtonText("Open File");

        chooseDirectory.showOpenDialog(null);
        currentFile = chooseDirectory.getSelectedFile().getAbsolutePath();
        if(currentFile != null)
            filenameTxt.setText(chooseDirectory.getSelectedFile().getName());
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
        if(source == numSkippedSpin)
        {
            numSkipped = Integer.valueOf(numSkippedSpin.getValue().toString()).intValue();
        }
        setValues();
    }

    public void skipHeader(StsAsciiFile asciiFile)
    {
        try
        {
            for (int i = 0; i < numSkipped; i++)
                asciiFile.readLine();
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING,  e.getMessage());
            return;
        }
    }

    public int getNumSkippedRows() { return numSkipped; }
    public boolean setValues()
    {
        if(currentFile == null)
            return false;

        StsFile file = StsFile.constructor(currentFile);
        StsAsciiFile asciiFile = new StsAsciiFile(file);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        try
        {
            fileTextArea.setText("");
            skipHeader(asciiFile);
            String line = asciiFile.readLine();
            int nLinesRead = 1;
            while((line != null) && (nLinesRead < 100))
            {
                appendLine(line);
                line = asciiFile.readLine();
                nLinesRead++;
            }
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "File header read error for " +
                    currentFile + ": " + e.getMessage());
            return false;
        }

        return true;
    }
    public void appendLine(String line)
    {
        fileTextArea.append(line + '\n');
    }
    public String getHeaderFile()
    {
        return currentFile;
    }
    public String getHeaderFilePath()
    {
        if(currentFile != null)
        {
            File file = new File(currentFile);
            return file.getPath();
        }
        else
        {
            return wizard.getModel().getProject().getRootDirString();
        }
    }
    public int getNumberRowsToSkip() { return numSkipped; }
}
