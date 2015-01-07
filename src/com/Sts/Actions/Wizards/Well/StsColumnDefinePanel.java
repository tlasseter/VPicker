package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.event.*;
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

public class StsColumnDefinePanel extends StsJPanel implements ActionListener, ChangeListener {
    private StsWellWizard wizard;
    private StsColumnDefine wizardStep;

    static final String[][][] labels = new String[StsWellWizard.NUM_FILE_TYPES][][];
    static final int[] numRequiredCols = new int[] {9, 3, 0, 2, 2};

    // Header File
    static final String[] uwiLabels = new String[] {"UWI","API","UWINO","APINO","UWI_NO","API_NO"};
    static final String[] nameLabels = new String[] {"NAME","WELLNAME","LABEL","LBL"};
    static final String[] xLabels = new String[] {"X", "XCOOR", "SURFX"};
    static final String[] yLabels = new String[] {"Y", "YCOOR", "SURFY"};
    static final String[] symbolLabels = new String[] {"SYMBOL","SYMB","SYM"};
    static final String[] kbLabels = new String[] {"KB", "KB_ELEV", "ELEV_KB","ELEVKB","KBELEV"};
    static final String[] grdLabels = new String[] {"GRD", "ELEV_GR", "GR_ELEV", "GRELEV", "ELV_GR", "ELEVGR","ELVGR","ELEV"};
    static final String[] tvdLabels = new String[] {"TD", "TVD", "DEPTH"};
    static final String[] datumLabels = new String[] {"DATUM", "DTM"};
    static final String[][] headerFileLabels = new String[9][];

    // Survey File
    static final String[] mdLabels = new String[] {"MD","MDEPTH","MDEP"};
    static final String[] azimLabels = new String[] {"AZIM", "AZIMUTH", "AZ", "AZM"};
    static final String[] dipLabels = new String[] {"DIP", "DP", "INCL", "INCLINATION"};
    static final String[][] surveyFileLabels = new String[3][];

    // Time-Depth File - uses MD from above
    static final String[] timeLabels = new String[] {"TIME", "TWT", "OWT"};
    static final String[][] tdFileLabels = new String[2][];

    // Tops File - uses MD from above
    static final String[] topLabels = new String[] {"MKR", "MARKER", "TOP", "TOPS", "MKRS", "MARKERS"};
    static final String[][] topFileLabels = new String[2][];

    StsGroupBox typeSelectionBox = new StsGroupBox();
    JRadioButton headerRadio = new JRadioButton();
    JLabel headerLabel = new JLabel("Header File");

    JRadioButton surveyRadio = new JRadioButton();
    JLabel surveyLabel = new JLabel("Survey File");

    JRadioButton topsRadio = new JRadioButton();
    JLabel topsLabel = new JLabel("Tops File");

    JRadioButton tdsRadio = new JRadioButton();
    JLabel tdsLabel = new JLabel("Time-Depth File");

    ButtonGroup fileTypeGroup = new ButtonGroup();

    JPanel jPanel1 = new JPanel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JScrollPane jScrollPane1 = new JScrollPane();
    JLabel fileViewLbl = new JLabel();
    JTextArea fileTextArea = new JTextArea();
    Font font = new Font("Serif", 1, 10);
    int height = 500;
    int width = 500;
    private StsGroupBox colDefPanel = new StsGroupBox();
    int numCols = 0;

    // Spinners for header file
    JLabel uwiLabel = new JLabel("UWI:");
    JLabel nameLabel = new JLabel("Name:");
    JLabel xLabel = new JLabel("Surface X:");
    JLabel yLabel = new JLabel("Surface Y:");
    JLabel symbolLabel = new JLabel("Symbol:");
    JLabel kbLabel = new JLabel("KB Elevation:");
    JLabel grdLabel = new JLabel("Elevation:");
    JLabel tdLabel = new JLabel("Total Depth:");
    JLabel datumLabel = new JLabel("Datum:");
    private JLabel[] headerSpinnerLabels = new JLabel[] {uwiLabel, nameLabel, xLabel,yLabel, symbolLabel, kbLabel,grdLabel,tdLabel,datumLabel};

    // Spinners for time-depth file
    JLabel depthLabel = new JLabel("Depth:");
    JLabel timeLabel = new JLabel("Time:");
    private JLabel[] tdSpinnerLabels = new JLabel[] {depthLabel,timeLabel};

    // Spinners for top file
    JLabel topLabel = new JLabel("Top:");
    JLabel topDepthLabel = new JLabel("Depth:");
    private JLabel[] topSpinnerLabels = new JLabel[] {topLabel,topDepthLabel};

    // Spinners for Survey File
    JLabel mdLabel = new JLabel("Depth:");
    JLabel azimuthLabel = new JLabel("Azimuth:");
    JLabel dipLabel = new JLabel("Dip:");
    private JLabel[] surveySpinnerLabels = new JLabel[] {mdLabel,azimuthLabel,dipLabel};

    private JSpinner[][] spinners = new JSpinner[StsWellWizard.NUM_FILE_TYPES][];
    private SpinnerModel[][] spinnerModels = new SpinnerModel[StsWellWizard.NUM_FILE_TYPES][];
    private JLabel[][] spinnerLabels = new JLabel[StsWellWizard.NUM_FILE_TYPES][];

    GridBagLayout gridBagLayout1 = new GridBagLayout();

    private String currentFile = null;

    public StsColumnDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsWellWizard)wizard;
        this.wizardStep = (StsColumnDefine)wizardStep;
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
        for(int i=0; i<StsWellWizard.NUM_FILE_TYPES; i++)
        {
            if(numRequiredCols[i] == 0) continue;
            spinners[i] = new JSpinner[numRequiredCols[i]];
            spinnerModels[i] = new SpinnerModel[numRequiredCols[i]];
            for(int j=0; j<numRequiredCols[i]; j++)
            {
                spinners[i][j] = new JSpinner();
                spinnerModels[i][j] = new SpinnerNumberModel(0, 0, 50, 1);
            }
        }
        spinnerLabels[StsWellWizard.SURVEY] = surveySpinnerLabels;
        spinnerLabels[StsWellWizard.TD] = tdSpinnerLabels;
        spinnerLabels[StsWellWizard.TOPS] = topSpinnerLabels;
        spinnerLabels[StsWellWizard.HEADER] = headerSpinnerLabels;

        headerFileLabels[StsWellKeywordIO.UWI] = uwiLabels;
        headerFileLabels[StsWellKeywordIO.NAME] = nameLabels;
        headerFileLabels[StsWellKeywordIO.MX] = yLabels;
        headerFileLabels[StsWellKeywordIO.MY] = xLabels;
        headerFileLabels[StsWellKeywordIO.SYMBOL] = symbolLabels;
        headerFileLabels[StsWellKeywordIO.KB] = kbLabels;
        headerFileLabels[StsWellKeywordIO.GRD] = grdLabels;
        headerFileLabels[StsWellKeywordIO.TVD] = tvdLabels;
        headerFileLabels[StsWellKeywordIO.DATUM] = datumLabels;

        surveyFileLabels[StsWellKeywordIO.MD] = mdLabels;
        surveyFileLabels[StsWellKeywordIO.AZIM] = azimLabels;
        surveyFileLabels[StsWellKeywordIO.DIP] = dipLabels;

        tdFileLabels[StsWellKeywordIO.TD_TIME] = timeLabels;
        tdFileLabels[StsWellKeywordIO.MD] = mdLabels;

        topFileLabels[StsWellKeywordIO.TOP] = topLabels;
        topFileLabels[StsWellKeywordIO.MD] = mdLabels;

        labels[StsWellWizard.SURVEY] = surveyFileLabels;
        labels[StsWellWizard.TOPS] = topFileLabels;
        labels[StsWellWizard.TD] = tdFileLabels;
        labels[StsWellWizard.HEADER] = headerFileLabels;

        initializeColOrdering();
        setValues();
    }

    private void initializeColOrdering()
    {
        // Examine Header
        String[] headerTokens = StsMultiWellImport.getHeaderRow(wizard.getFilename(StsWellWizard.HEADER), wizard.getNumberRowsToSkip(StsWellWizard.HEADER));
        determineDefaultCols(headerTokens, StsWellWizard.HEADER);

        // Examine Survey
        headerTokens = StsMultiWellImport.getHeaderRow(wizard.getFilename(StsWellWizard.SURVEY), wizard.getNumberRowsToSkip(StsWellWizard.SURVEY));
        determineDefaultCols(headerTokens, StsWellWizard.SURVEY);

        // Examine Time-Depth
        headerTokens = StsMultiWellImport.getHeaderRow(wizard.getFilename(StsWellWizard.TD), wizard.getNumberRowsToSkip(StsWellWizard.TD));
        determineDefaultCols(headerTokens, StsWellWizard.TD);

        // Examine Tops
        headerTokens = StsMultiWellImport.getHeaderRow(wizard.getFilename(StsWellWizard.TOPS), wizard.getNumberRowsToSkip(StsWellWizard.TOPS));
        determineDefaultCols(headerTokens, StsWellWizard.TOPS);
    }

    public void determineDefaultCols(String[] hdrTokens, byte type)
    {
        boolean foundIt = false;
        int loc = 0;
        if(hdrTokens == null) return;
        for(int i=0; i<labels[type].length; i++)
        {
            foundIt = false;
            loc = 0;
            for(int k=0; k<labels[type][i].length; k++)
            {
                String label = labels[type][i][k];
                for(int j=0; j<hdrTokens.length; j++)
                {
                    if(label.indexOf(hdrTokens[j].toUpperCase()) >= 0)
                    {
                        foundIt = true;
                        loc = j;
                        break;
                    }
                }
                if(foundIt)
                {
                    spinners[type][i].setValue(new Integer(loc));
                    break;
                }
            }

        }
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout1);
        jPanel1.setLayout(gridBagLayout2);

        headerRadio.addActionListener(this);
        typeSelectionBox.addToRow(headerRadio,1,0.0);
        typeSelectionBox.addEndRow(headerLabel,1,0.0);

        surveyRadio.addActionListener(this);
        typeSelectionBox.addToRow(surveyRadio,1,0.0);
        typeSelectionBox.addEndRow(surveyLabel,1,0.0);

        tdsRadio.addActionListener(this);
        typeSelectionBox.addToRow(tdsRadio,1,0.0);
        typeSelectionBox.addEndRow(tdsLabel,1,0.0);

        topsRadio.addActionListener(this);
        typeSelectionBox.addToRow(topsRadio,1,0.0);
        typeSelectionBox.addEndRow(topsLabel,1,0.0);

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

        for(int i=0; i<spinners.length; i++)
        {
            if(spinners[i] == null) continue;
            for(int j=0; j<spinners[i].length; j++)
            {
                spinners[i][j].addChangeListener(this);
                spinners[i][j].setModel(spinnerModels[i][j]);
            }
        }
        headerRadio.setSelected(true);
        configureColPanel(StsWellWizard.SURVEY);

        gbc.fill = gbc.HORIZONTAL;
        gbc.gridheight = 1;
        gbc.weighty = 0.0;
        gbc.anchor = gbc.NORTH;
        this.addEndRow(typeSelectionBox,3,1.0);

        gbc.fill = gbc.BOTH;
        gbc.gridheight = 5;
        gbc.weighty = 1.0;
        this.addToRow(colDefPanel,1,0.1);
        this.addEndRow(jPanel1,2,1.0);

        fileTypeGroup.add(surveyRadio);
        fileTypeGroup.add(headerRadio);
        fileTypeGroup.add(topsRadio);
        fileTypeGroup.add(tdsRadio);
    }

    public void configureColPanel(byte type)
    {
        colDefPanel.removeAll();
        if(spinners[type] == null) return;
        colDefPanel.gbc.anchor = gbc.NORTH;
        for(int i=0; i<spinners[type].length; i++)
        {
            colDefPanel.addToRow(spinnerLabels[type][i]);
            colDefPanel.addEndRow(spinners[type][i]);
        }
        wizard.rebuild();
        this.repaint();
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == surveyRadio)
            configureColPanel(StsWellWizard.SURVEY);
        else if(source == headerRadio)
            configureColPanel(StsWellWizard.HEADER);
        else if(source == topsRadio)
            configureColPanel(StsWellWizard.TOPS);
        else if(source == tdsRadio)
            configureColPanel(StsWellWizard.TD);

        currentFile = wizard.getFilename(getSelectedFileType());
        setValues();
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
        int type = getSelectedFileType();

        // ToDo: Change the column lable or highlighting to show user change.
    }

    private byte getSelectedFileType()
    {
        if(surveyRadio.isSelected()) return StsWellWizard.SURVEY;
        else if(topsRadio.isSelected()) return StsWellWizard.TOPS;
        else if(headerRadio.isSelected()) return StsWellWizard.HEADER;
        else if(tdsRadio.isSelected()) return StsWellWizard.TD;
        else
        {
            headerRadio.setSelected(true);
            return StsWellWizard.HEADER;
        }
    }

    public void skipHeader(StsAsciiFile asciiFile)
    {
        try
        {
            for (int i = 0; i < wizard.getNumberRowsToSkip(getSelectedFileType()); i++)
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

    public int getNumRowsToSkip(int type)
    {
        return wizard.getNumberRowsToSkip(type);
    }

    public int[] getColOrder(byte type)
    {
        int[] colOrder = new int[spinners[type].length];
        int numCols = spinners[type].length;

        for(int i=0; i<spinners[StsWellWizard.HEADER].length; i++)
            colOrder[i] = Integer.parseInt((String)spinnerModels[StsWellWizard.HEADER][i].getValue());

        return colOrder;
    }
}
