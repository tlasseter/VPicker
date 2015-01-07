package com.Sts.UI;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import javax.swing.table.*;
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

public class StsTextHeaderEditDialog extends JDialog implements ActionListener, ChangeListener
{
    private ButtonGroup btnGroup1 = new ButtonGroup();
    public StsModel model = null;
    Font defaultFont = new Font("Dialog",0,11);
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    TitledBorder titledBorder1;
    SpinnerListModel spinnerModel = null;
    String[] list = null;
    String header = null;
    private JTable headerTable = null;
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JButton okBtn = new JButton();
    JButton cancelBtn = new JButton();

    static byte OK = 0;
    static byte CANCELED = 1;
    private byte mode = OK;

    public StsTextHeaderEditDialog(Frame frame, String header, boolean modal)
    {
        super(frame,"SegY Export Text Header Viewer / Editor", modal);
        this.header = header;
        //this.setLocationRelativeTo(frame);
		StsToolkit.centerComponentOnScreen(this);

        try
        {
            initialize();
            jbInit();
            pack();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        titledBorder1 = new TitledBorder("");
        this.setModal(true);
        this.getContentPane().setLayout(gridBagLayout2);

        JScrollPane scrollPane = new JScrollPane(headerTable);
        okBtn.setText("Ok");
        okBtn.addActionListener(this);
        cancelBtn.setText("Cancel");
        cancelBtn.addActionListener(this);
        headerTable.setSelectionBackground(SystemColor.inactiveCaption);
        headerTable.setSelectionForeground(Color.black);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setPreferredSize(new Dimension(330, 670));
        scrollPane.setRequestFocusEnabled(true);
        scrollPane.getViewport().add(headerTable);
        headerTable.setCellSelectionEnabled(true);
        headerTable.setDragEnabled(true);
//        headerTable.setSelectionBackground();

        this.getContentPane().add(scrollPane,           new GridBagConstraints(0, 0, 2, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 5, 0, 5), 0, 0));
        this.getContentPane().add(okBtn,            new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 5, 5, 0), 30, 0));
        this.getContentPane().add(cancelBtn,          new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 5), 30, 0));
    }

    private void initialize()
    {
        final String[] columnNames = {" Row ","Comment (cols 5-80)",};
        final Object[][] rowData = new Object[40][2];

        for(int i=0; i<40; i++)
        {
            rowData[i][0] = header.substring(i*80, (i*80)+4);
            rowData[i][1] = header.substring((i*80)+4, (i*80)+80);
        }

        TableModel myModel = new AbstractTableModel() {
            public String getColumnName(int col) {
                return columnNames[col].toString();
            }
            public int getRowCount() { return rowData.length; }
            public int getColumnCount() { return columnNames.length; }
            public Object getValueAt(int row, int col) {
                return rowData[row][col];
            }
            public boolean isCellEditable(int row, int col)
                { return true; }
            public void setValueAt(Object value, int row, int col) {
                rowData[row][col] = value;
                fireTableCellUpdated(row, col);
            }
        };
        headerTable = new JTable(myModel);

        headerTable.getColumnModel().getColumn(0).setMinWidth(50);
        headerTable.getColumnModel().getColumn(0).setMaxWidth(50);
        headerTable.getColumnModel().getColumn(1).setWidth(100);
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == okBtn)
        {
            mode = OK;
            String line = null;

            for(int i=0; i<40; i++)
            {
                line = ((String)headerTable.getValueAt(i,0)).substring(0,4) + (String)headerTable.getValueAt(i,1);
                line = StsStringUtils.padClipString(line, 80);
                if(i == 0)
                    header = line;
                else
                    header = header + line;
            }
            setVisible(false);
        }
        else if(source == cancelBtn)
        {
            mode = CANCELED;
            setVisible(false);
        }
    }

    public byte getMode() { return mode;  }

    public String getHeaderText()
    {
        return header;
    }
}
