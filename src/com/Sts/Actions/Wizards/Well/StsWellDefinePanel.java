package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.text.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsWellDefinePanel extends JPanel implements ActionListener, ChangeListener
{
        private ButtonGroup zGroup = new ButtonGroup();
        private JLabel jLabel1 = new JLabel();
        DefaultTableModel segyTableModel = new DefaultTableModel();
        private StsTablePanel fileTable = new StsTablePanel();
        private JLabel jLabel4 = new JLabel();
        private JLabel jLabel5 = new JLabel();
        private JLabel jLabel6 = new JLabel();
        private JComboBox fileCombo = new JComboBox();
        private String currentFile = null;
        private float currentRowNum = -1;
        private int headerFormat = StsSEGYFormat.INT4;
        private DecimalFormat labelFormat = new DecimalFormat("####.#");

        private StsWizard wizard;
        private StsDefineUtWells wizardStep;

        private StsProperties textHeaderFormatTable;

        int selectedCol, xCol, yCol;
        int numCols = 0;
        double zValue = 0.0f;

        JPanel colDefPanel = new JPanel();
        JLabel jLabel2 = new JLabel();
        JSpinner xSpin = new JSpinner();
        private SpinnerModel xSpinModel = null;
        JSpinner ySpin = new JSpinner();
        private SpinnerModel ySpinModel = null;

        JLabel jLabel8 = new JLabel();
        JLabel jLabel9 = new JLabel();
        JLabel jLabel10 = new JLabel();
        JLabel jLabel11 = new JLabel();
        JLabel jLabel12 = new JLabel();
        JPanel jPanel1 = new JPanel();
        JPanel jPanel2 = new JPanel();
        JButton okBtn = new JButton();
        JButton cancelBtn = new JButton();
        GridBagLayout gridBagLayout3 = new GridBagLayout();
        JLabel jLabel13 = new JLabel();
        StsFloatFieldBean zBean = new StsFloatFieldBean();
        GridBagLayout gridBagLayout1 = new GridBagLayout();
        GridBagLayout gridBagLayout2 = new GridBagLayout();
        GridBagLayout gridBagLayout4 = new GridBagLayout();

        public StsWellDefinePanel(StsWizard wizard, StsDefineUtWells wizardStep)
        {
            this.wizard = wizard;
            this.wizardStep = wizardStep;
            try
            {
                jbInit();
            }
            catch(Exception ex) {
                ex.printStackTrace();
            }
        }

        public void initialize()
        {
            xCol = 3;
            yCol = 4;
            setColumnHeadings();
        }

        void jbInit() throws Exception
        {
            this.setLayout(gridBagLayout4);

            xSpinModel = new SpinnerNumberModel(4, 1, 20, 1);
            ySpinModel = new SpinnerNumberModel(5, 1, 20, 1);
            zBean.initialize(this, "Z", true, "Kelly Bushing Elevation: ");

//            fileTable.setPreferredSize(new Dimension(500, 500));
            fileTable.setTitle("Well File:");

            colDefPanel.setBorder(BorderFactory.createEtchedBorder());
            colDefPanel.setLayout(gridBagLayout1);
            jLabel8.setHorizontalAlignment(SwingConstants.RIGHT);
            jLabel8.setText("X:");
            jLabel9.setHorizontalAlignment(SwingConstants.RIGHT);
            jLabel9.setText("Y:");
            xSpin.addChangeListener(this);
            xSpin.setModel(xSpinModel);
            ySpin.addChangeListener(this);
            ySpin.setModel(ySpinModel);
            jPanel1.setBorder(BorderFactory.createEtchedBorder());
            jPanel1.setLayout(gridBagLayout2);
            jPanel2.setBorder(BorderFactory.createEtchedBorder());
            jPanel2.setLayout(gridBagLayout3);
            okBtn.setText("Ok");
            okBtn.addActionListener(this);
            cancelBtn.setText("Cancel");
            cancelBtn.addActionListener(this);
            jLabel13.setText("Kelly Bushing Elevation:");
            jLabel13.setHorizontalAlignment(SwingConstants.RIGHT);
            zBean.setText("0.0");
        jPanel2.add(cancelBtn, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 236, 5, 4), 8, 0));
        jPanel2.add(okBtn, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 5, 0), 28, 0));
        colDefPanel.add(jLabel12,  new GridBagConstraints(0, 0, GridBagConstraints.REMAINDER, GridBagConstraints.REMAINDER, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
            colDefPanel.add(jLabel8,   new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 6, 2, 0), 0, 0));
            colDefPanel.add(xSpin,   new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 0), 0, 0));
            colDefPanel.add(jLabel13,    new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 25, 2, 0), 0, 0));
            colDefPanel.add(ySpin,   new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(2, 2, 2, 0), 0, 0));
            colDefPanel.add(jLabel9,   new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(9, 25, 3, 0), 0, 0));
        colDefPanel.add(zBean,  new GridBagConstraints(5, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 3, 2), 4, 5));
            jPanel1.add(fileTable,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 5, 0, 0), -117, -329));
            jPanel1.add(colDefPanel,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 0, 0), 1, 0));
            this.add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 2, 0, 6), 0, 0));
        this.add(jPanel2,   new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 2, 5, 6), 0, 1));

        if(!setValues())
            return;
        }

        public void actionPerformed(ActionEvent e)
        {
            Object source = e.getSource();
            if(source == okBtn)
            {
                wizardStep.panel.exitFileView();
            }
            else if(source == cancelBtn)
            {
                xCol = -1;
                yCol = -1;
                wizardStep.panel.exitFileView();
            }
        }

        public void stateChanged(ChangeEvent e)
        {
            Object source = e.getSource();
            if(source == xSpin)
            {
                xCol = Integer.valueOf(xSpin.getValue().toString()).intValue() - 1;
            }
            else if(source == ySpin)
            {
                yCol = Integer.valueOf(ySpin.getValue().toString()).intValue() - 1;
            }
            setColumnHeadings();
        }

        public boolean setValues()
        {
            String[] tokens = null;
            int nLines = 0;

            if(wizardStep.panel.getSelectedFile() == null)
                return false;

            String dir = StsWellImport.getCurrentDirectory() + File.separator;
            String filename = wizardStep.panel.getSelectedFile();
            StsFile file = StsFile.constructor(dir, filename);
            StsAsciiFile asciiFile = new StsAsciiFile(file);
            if(!asciiFile.openReadWithErrorMessage()) return false;
         
            double value;
            try
            {
                while(true)
                {
                    tokens = asciiFile.getTokens();
                    try { value = Double.parseDouble(tokens[0]); break;}
                    catch(Exception e) { continue; }
                }

                int maxTokens = tokens.length;
                Object[] row = new Object[maxTokens];
                numCols = tokens.length;
                if (numCols < 5)
                {
                    new StsMessage(wizard.getModel().win3d, StsMessage.ERROR, "Insufficient entries for line: " + asciiFile.getLine() + " need 5 columns.");
                    xCol = -1;
                    yCol = -1;
                    wizardStep.panel.exitFileView();
                    return false;
                }
                setColumnHeadings();
                fileTable.removeAllRows();

                int maxLineErrors = 10;
                int nLineErrors = 0;
                while(tokens != null)
                {
                    if(++nLines > 100)
                        break;
                    int nTokens = tokens.length;
                    if (nTokens < 5)
                    {
                        if(nLineErrors < maxLineErrors)
                        {
                            nLineErrors++;
                            if (nLineErrors < maxLineErrors) StsMessageFiles.errorMessage("Insufficient entries for line " + nLines + ": " + asciiFile.getLine());
                            else if(nLineErrors == maxLineErrors) StsMessageFiles.errorMessage("More errors...insufficient entries.");
                        }
                        continue;
                    }
                    if(tokens.length > maxTokens && nLineErrors < maxLineErrors)
                    {
                        nLineErrors++;
                        if(nLineErrors < maxLineErrors) StsMessageFiles.errorMessage("Number of tokens exceeds " + maxTokens + " line " + nLines + ": " + asciiFile.getLine());
                        else if(nLineErrors == maxLineErrors) StsMessageFiles.errorMessage("More errors...lines have more than " + maxTokens + " tokens.");
                    }
                    int nTokensToLoad = Math.min(tokens.length, maxTokens);
                    for(int j = 0; j < nTokensToLoad; j++)
                    {
                        try
                        {
                            value = Double.parseDouble(tokens[j]);
                            row[j] = labelFormat.format(value);
                        }
                        catch(Exception e)
                        {
                            if(j == 0)
                                break;
                            else
                                row[j] = tokens[j];
                        }
                    }
                    fileTable.addRow(row);
                    tokens = asciiFile.getTokens();
                }
                return true;
            }
            catch(IOException e)
            {
                new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "File read error for " +
                        filename + " line " + nLines + ": " + asciiFile.getLine());
                return false;
            }
        }

        public void setColumnHeadings()
        {
            Object[] colNames = new Object[numCols];
            for(int i=0; i<numCols; i++)
                colNames[i] = "Unknown";

            if(xCol > numCols-1)
                xCol = numCols-1;
            if(yCol > numCols-1)
                yCol = numCols-1;

            colNames[0] = "MD";
            colNames[1] = "Drift";
            colNames[2] = "Azimuth";
            colNames[xCol] = "X";
            colNames[yCol] = "Y";

            fileTable.changeColumnHeadings(colNames);
        }

        public int getYCol()
        {
            return yCol;
        }

        public int getXCol()
        {
            return xCol;
        }

        public double getZ()
        {
            return zValue;
        }
        public void setZ(double z)
        {
            zValue = z;
        }

}
