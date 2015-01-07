package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.UI.Beans.*;
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

public class StsOtherFilesPanel extends StsJPanel implements ActionListener, ChangeListener {
    private StsWellWizard wizard;
    private StsOtherFiles wizardStep;

    private JFileChooser chooseDirectory = null;

    StsGroupBox surveyFileSelectionBox = new StsGroupBox();
    JRadioButton surveyRadio = new JRadioButton();
    JLabel surveyLabel = new JLabel("Survey File:");
    JButton surveyBtn = new JButton();
    JTextField surveyFileTxt = new JTextField();
    int numSurveyRowsSkipped = 0;

    StsGroupBox topsFileSelectionBox = new StsGroupBox();
    JRadioButton topsRadio = new JRadioButton();
    JLabel topsLabel = new JLabel("Tops File:");
    JButton topsBtn = new JButton();
    JTextField topsFileTxt = new JTextField();
    int numTopsRowsSkipped = 0;

    private String[] logFiles = null;
    StsGroupBox logsFileSelectionBox = new StsGroupBox();
    JRadioButton logsRadio = new JRadioButton();
    JLabel logsLabel = new JLabel("Log File:");
    JButton logsBtn = new JButton();
    StsComboBoxFieldBean logsFileBean = new StsComboBoxFieldBean(this, "logFile", "Log Files:", logFiles);
    int numLogsRowsSkipped = 0;

    StsGroupBox tdsFileSelectionBox = new StsGroupBox();
    JRadioButton tdsRadio = new JRadioButton();
    JLabel tdsLabel = new JLabel("Time-Depth File:");
    JButton tdsBtn = new JButton();
    JTextField tdsFileTxt = new JTextField();
    int numTdsRowsSkipped = 0;

    ButtonGroup fileTypeGroup = new ButtonGroup();

    JPanel jPanel1 = new JPanel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JScrollPane jScrollPane1 = new JScrollPane();
    JLabel fileViewLbl = new JLabel();
    JTextArea fileTextArea = new JTextArea();
    Font font = new Font("Serif", 1, 10);
    int height = 500;
    int width = 500;
    private StsGroupBox colDefPanel = new StsGroupBox("# Header Rows");
    int numCols = 0;
    //int numSkipped = 0;
    JLabel jLabel2 = new JLabel();
    JSpinner numSkippedSpin = new JSpinner();
    private SpinnerModel numSkippedSpinModel = null;

    GridBagLayout gridBagLayout1 = new GridBagLayout();

    private String currentFile = null;
    private String currentSurveyFile = null;
    private String currentTopsFile = null;
    private String currentLogsFile = null;
    private String currentTdFile = null;

    private byte selectedFileType = StsWellWizard.SURVEY;

    public StsOtherFilesPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsWellWizard)wizard;
        this.wizardStep = (StsOtherFiles)wizardStep;
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

        surveyBtn.setMargin(new Insets(0, 2, 0, 2));
        surveyBtn.setText("Select...");
        surveyBtn.addActionListener(this);
        surveyRadio.addActionListener(this);
        surveyFileTxt.setText("Survey File");
        surveyFileSelectionBox.addToRow(surveyRadio,1,0.0);
        surveyFileSelectionBox.addToRow(surveyLabel,1,0.0);
        surveyFileSelectionBox.gbc.fill = gbc.HORIZONTAL;
        surveyFileSelectionBox.addToRow(surveyFileTxt,2,1.0);
        surveyFileSelectionBox.gbc.fill = gbc.NONE;
        surveyFileSelectionBox.addEndRow(surveyBtn,1,0.0);

        topsBtn.setMargin(new Insets(0, 2, 0, 2));
        topsBtn.setText("Select...");
        topsBtn.addActionListener(this);
        topsRadio.addActionListener(this);
        topsFileTxt.setText("Tops File");
        topsFileSelectionBox.addToRow(topsRadio,1,0.0);
        topsFileSelectionBox.addToRow(topsLabel,1,0.0);
        topsFileSelectionBox.gbc.fill = gbc.HORIZONTAL;
        topsFileSelectionBox.addToRow(topsFileTxt,2,1.0);
        topsFileSelectionBox.gbc.fill = gbc.NONE;
        topsFileSelectionBox.addEndRow(topsBtn,1,0.0);

        logsBtn.setMargin(new Insets(0, 2, 0, 2));
        logsBtn.setText("Select...");
        logsBtn.addActionListener(this);
        logsRadio.addActionListener(this);
        logsFileSelectionBox.addToRow(logsRadio,1,0.0);
        logsFileSelectionBox.gbc.fill = gbc.HORIZONTAL;
        logsFileSelectionBox.addToRow(logsFileBean,3,1.0);
        logsFileSelectionBox.gbc.fill = gbc.NONE;
        logsFileSelectionBox.addEndRow(logsBtn,1,0.0);

        tdsBtn.setMargin(new Insets(0, 2, 0, 2));
        tdsBtn.setText("Select...");
        tdsBtn.addActionListener(this);
        tdsRadio.addActionListener(this);
        tdsFileTxt.setText("Time-Depth File");
        tdsFileSelectionBox.addToRow(tdsRadio,1,0.0);
        tdsFileSelectionBox.addToRow(tdsLabel,1,0.0);
        tdsFileSelectionBox.gbc.fill = gbc.HORIZONTAL;
        tdsFileSelectionBox.addToRow(tdsFileTxt,2,1.0);
        tdsFileSelectionBox.gbc.fill = gbc.NONE;
        tdsFileSelectionBox.addEndRow(tdsBtn,1,0.0);

        fileTextArea.setText("First 100 lines of file will be isVisible here on selection...");
        fileViewLbl.setFont(new Font("Serif", 1, 12));
        fileViewLbl.setHorizontalAlignment(SwingConstants.CENTER);
        fileViewLbl.setText("File View");
        fileTextArea.setBackground(Color.lightGray);
        fileTextArea.setEditable(false);
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
        this.addEndRow(surveyFileSelectionBox,3,1.0);
        this.addEndRow(logsFileSelectionBox,3,1.0);
        this.addEndRow(topsFileSelectionBox,3,1.0);
        this.addEndRow(tdsFileSelectionBox,3,1.0);

        gbc.fill = gbc.BOTH;
        gbc.gridheight = 10;
        gbc.weighty = 1.0;
        this.addToRow(colDefPanel,1,0.0);
        this.addEndRow(jPanel1,2,1.0);

        surveyRadio.setSelected(true);

        fileTypeGroup.add(surveyRadio);
        fileTypeGroup.add(logsRadio);
        fileTypeGroup.add(topsRadio);
        fileTypeGroup.add(tdsRadio);
    }

    public String getLogFile()
    {
        return currentLogsFile;
    }

    public void setLogFile(String logfile)
    {
        for(int i=0; i< logFiles.length; i++)
        {
//            if(logFiles[i].contains(logfile))
            {
                currentFile = logFiles[i];
                break;
            }
        }
        currentLogsFile = currentFile;
        logsRadio.setSelected(true);
        numSkippedSpin.setEnabled(false);
        setValues();
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        int numRowsSkipped = 0;
        try
        {
            if(source == surveyBtn)
            {
                surveyRadio.setSelected(true);
                directoryBrowse();
                surveyFileTxt.setText(StsFile.getFilenameFromPathname(currentFile));
                currentSurveyFile = currentFile;
                numRowsSkipped =  numSurveyRowsSkipped;
                numSkippedSpin.setEnabled(true);
            }

            else if(source == topsBtn)
            {
                topsRadio.setSelected(true);
                directoryBrowse();
                topsFileTxt.setText(StsFile.getFilenameFromPathname(currentFile));
                currentTopsFile = currentFile;
                numRowsSkipped =  numTopsRowsSkipped;
                numSkippedSpin.setEnabled(true);
            }

            else if(source == logsBtn)
            {
                logsRadio.setSelected(true);
                directoryBrowse();
                logsFileBean.setSelectedIndex(0);
                currentLogsFile = currentFile;
                numRowsSkipped =  numLogsRowsSkipped;
                numSkippedSpin.setEnabled(false);
            }

            else if(source == tdsBtn)
            {
                tdsRadio.setSelected(true);
                directoryBrowse();
                tdsFileTxt.setText(StsFile.getFilenameFromPathname(currentFile));
                currentTdFile = currentFile;
                numRowsSkipped =  numTdsRowsSkipped;
                numSkippedSpin.setEnabled(true);
            }

            else if(source == surveyRadio)
            {
                numRowsSkipped =  numSurveyRowsSkipped;
                currentFile = currentSurveyFile;
            }

            else if(source == logsRadio)
            {
                numRowsSkipped =  numLogsRowsSkipped;
                currentFile = currentLogsFile;
            }

            else if(source == topsRadio)
            {
                numRowsSkipped =  numTopsRowsSkipped;
                currentFile = currentTopsFile;
            }

            else if(source == tdsRadio)
            {
                numRowsSkipped =  numTdsRowsSkipped;
                currentFile = currentTdFile;
            }
        }
        catch(Exception ex)
        {
            StsException.outputException(ex, StsException.WARNING);
        }
        numSkippedSpinModel.setValue(new Integer(numRowsSkipped));
        if(currentFile != null)
            setValues();
    }

    public void directoryBrowse()
    {
        String currentDirectory = wizard.getDefaultFilePath();
        chooseDirectory = new JFileChooser(currentDirectory);

        try
        {
            if(getSelectedFileType() != StsWellWizard.LOGS)
            {
                chooseDirectory.setDialogTitle("Select or Enter Desired File and Press Open");
                chooseDirectory.setApproveButtonText("Open File");
                chooseDirectory.setMultiSelectionEnabled(false);
            }
            else
            {
                chooseDirectory.setDialogTitle("Multi-Select Desired Files and Press Open");
                chooseDirectory.setApproveButtonText("Open Files");
                chooseDirectory.setMultiSelectionEnabled(true);
            }
            chooseDirectory.setFileSelectionMode(JFileChooser.FILES_ONLY);
            chooseDirectory.showOpenDialog(null);

            if(getSelectedFileType() == StsWellWizard.LOGS)
            {
                logsFileBean.removeAll();
                File[] files = chooseDirectory.getSelectedFiles();
                logFiles = new String[files.length];
                for(int i=0; i<files.length; i++)
                {
                    logFiles[i] = files[i].getAbsolutePath();
                    logsFileBean.addItem(StsFile.getFilenameFromPathname(logFiles[i]));
                }
                currentFile = logFiles[0];
            }
            else
                currentFile = chooseDirectory.getSelectedFile().getAbsolutePath();
        }
        catch(Exception e)
        {
            StsException.outputException(e,StsException.WARNING);
        }
    }

    private byte getSelectedFileType()
    {
        if(surveyRadio.isSelected()) return StsWellWizard.SURVEY;
        else if(topsRadio.isSelected()) return StsWellWizard.TOPS;
        else if(logsRadio.isSelected()) return StsWellWizard.LOGS;
        else if(tdsRadio.isSelected()) return StsWellWizard.TD;
        else
        {
            surveyRadio.setSelected(true);
            return StsWellWizard.SURVEY;
        }
    }

    public int getNumRowsSkipped()
    {
        switch(getSelectedFileType())
        {
            case StsWellWizard.SURVEY:
                return numSurveyRowsSkipped;
            case StsWellWizard.TOPS:
                return numTopsRowsSkipped;
            case StsWellWizard.LOGS:
                return numLogsRowsSkipped;
            case StsWellWizard.TD:
                return numTdsRowsSkipped;
            default:
                surveyRadio.setSelected(true);
                return numSurveyRowsSkipped;
        }
    }
    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
        if(source == numSkippedSpin)
        {
            int numRowsSkipped = Integer.valueOf(numSkippedSpin.getValue().toString()).intValue();
            switch(getSelectedFileType())
            {
                case StsWellWizard.SURVEY:
                    numSurveyRowsSkipped = numRowsSkipped;
                    break;
                case StsWellWizard.TOPS:
                    numTopsRowsSkipped = numRowsSkipped;
                    break;
                case StsWellWizard.LOGS:
                    numLogsRowsSkipped = numRowsSkipped;
                    break;
                case StsWellWizard.TD:
                    numTdsRowsSkipped = numRowsSkipped;
                    break;
                default:
                    surveyRadio.setSelected(true);

            }
        }
        setValues();
    }

    public void skipHeader(StsAsciiFile asciiFile)
    {
        try
        {
            for (int i = 0; i < getNumRowsSkipped(); i++)
                asciiFile.readLine();
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING,  e.getMessage());
            return;
        }
    }

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
            while((line != null) && (nLinesRead < 25))
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
     public String getSurveyFile()
    {
        return currentSurveyFile;
    }
    public int getNumRowsToSkip(int type)
    {
        switch(type)
        {
            case StsWellWizard.LOGS:
                return 0;
            case StsWellWizard.SURVEY:
                return this.numSurveyRowsSkipped;
            case StsWellWizard.TOPS:
                return numTopsRowsSkipped;
            case StsWellWizard.TD:
                return numTdsRowsSkipped;
        }
        return 0;
    }
    public String[] getLogFiles()
    {
        return logFiles;
    }
     public String getTopsFile()
    {
        return currentTopsFile;
    }
     public String getTdFile()
    {
        return currentTdFile;
    }
}
