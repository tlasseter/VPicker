package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import javax.swing.event.TableModelEvent;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import java.lang.Exception;


/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsWellEditDialog extends JDialog implements ActionListener, ChangeListener, TableModelListener
{
    private ButtonGroup btnGroup1 = new ButtonGroup();
    public StsModel model = null;
    public StsWell well;
	transient StsObjectRefList logs;
	transient String[] logNames;
	transient StsWellViewModel wellModel=null;
	transient StsEditTd editTd = null;

    Font defaultFont = new Font("Dialog",0,11);
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    DateFormat dateFormat = DateFormat.getDateTimeInstance();
    private DecimalFormat labelFormat = new DecimalFormat("#####.##");
    private JScrollPane tablePane = new JScrollPane();
    private JScrollPane normalPane = new JScrollPane();
    private JScrollPane xplotPane = new JScrollPane();

    TitledBorder titledBorder1;
    SpinnerListModel spinnerModel = null;
    String[] list = null;
    private StsTablePanel logTable = new StsTablePanel();
    JTabbedPane tabPane = new JTabbedPane();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
  JLabel jLabel1 = new JLabel();
  JLabel jLabel2 = new JLabel();

    public StsWellEditDialog(StsWellViewModel wellModel, StsWell well, boolean modal)
    {
        super(wellModel.getParentFrame(),"Well & Log Viewer / Editor (" + well.getName() + ")", modal);
        this.well = well;
		this.wellModel = wellModel;
        this.setLocationRelativeTo(wellModel.getWellWindowPanel());
        this.setSize(new Dimension(250,400));
		editTd =  new StsEditTd(well, false);
		logs = well.getLogCurves();
		logNames = new String[logs.getSize()+1];
		StsLogCurve logCurve = (StsLogCurve)logs.getElement(0);

		StsPoint[] pts = well.getLineVertexPoints();
		for (int i =0; i < pts.length; i++) {
//			System.out.println("add point "+i+" "+pts[i]);
			editTd.addPoint(pts[i]);
		}
//        this.setLocationRelativeTo(wellModel.getContentPane());
//        this.setSize(new Dimension(400,600));

        try
        {
            jbInit();
            initialize(wellModel);
            pack();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        tablePane.setAutoscrolls(true);
        tablePane.setPreferredSize(new Dimension(700, 700));
        titledBorder1 = new TitledBorder("");
        this.setModal(false);
        this.getContentPane().setLayout(gridBagLayout2);

        tabPane.setTabPlacement(JTabbedPane.BOTTOM);
        tabPane.setBorder(BorderFactory.createEtchedBorder());
        jLabel1.setRequestFocusEnabled(true);
        jLabel1.setHorizontalAlignment(SwingConstants.CENTER);
        //jLabel1.setText("Add Freq Plot");
        jLabel2.setHorizontalAlignment(SwingConstants.CENTER);
        //jLabel2.setText("Add Log Crossplot");
        this.getContentPane().add(tabPane,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 7, 7, 8), 0, 0));
        tabPane.add(tablePane, "Table");
        //tabPane.add(normalPane, "Normalize");
        //normalPane.getViewport().add(jLabel1, null);
        //tabPane.add(xplotPane, "Crossplot");
        //xplotPane.getViewport().add(jLabel2, null);
        tablePane.getViewport().add(logTable, null);

        tabPane.setSelectedIndex(0);
    }

    private void initialize(StsWellViewModel wellModel)
    {
        int i, j;
		logs = well.getLogCurves();
        logNames = new String[logs.getSize()+1];
        Object[] row = new Object[((StsLogCurve)logs.getElement(0)).getValuesFloatVector().getSize()+1];

        logTable.setTitle("Well: " + well.getName());
        //logNames[0] = new String("MD");
        for (i = 0; i < logs.getSize(); i++)
            logNames[i] = ((StsLogCurve)logs.getElement(i)).getName();
        logTable.addColumns(logNames);

        for(j=0; j< ((StsLogCurve)logs.getElement(0)).getValuesFloatVector().getSize(); j++)
        {
            //row[0] = labelFormat.format(new Float(((StsLogCurve)logs.getElement(0)).getMDepthFloatVector().getElement(j)));
            for (i = 0; i < logs.getSize(); i++)
            {
                if(((StsLogCurve)logs.getElement(i)).getValuesFloatVector().getElement(j) != StsParameters.nullValue)
                    row[i] = labelFormat.format(new Float(((StsLogCurve)logs.getElement(i)).getValuesFloatVector().getElement(j)));
                else
                    row[i] = new String("-");
            }
            logTable.addRow(row);
        }
		logTable.table.getModel().addTableModelListener(this);
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
		System.out.println("performed "+e);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
		System.out.println("performed "+e);
    }
	private void saveIt(StsLogCurve curve)
	{
	String binaryDataDir = curve.getStsDirectory() + wellModel.getModel().getProject().getBinaryDirString();
	curve.getValueVector().checkWriteBinaryFile(binaryDataDir, false, true);
	}
	public void tableChanged(TableModelEvent e)
	{
		//ystem.out.println("chamged "+e);
		int col = e.getColumn();
		int rowb = e.getFirstRow();
		int rowe = e.getLastRow();
		for (int row=rowb; row <= rowe; row++)
		{
			float origval = ((StsLogCurve)logs.getElement(col)).getValuesFloatVector().getElement(row);
			float newval = origval;
			String s = (String) logTable.table.getValueAt(row, col);
			try {
				newval = Float.valueOf(s.trim()).floatValue();
				((StsLogCurve)logs.getElement(col)).getValuesFloatVector().setElement(row, newval);
				wellModel.display();
				saveIt ((StsLogCurve)logs.getElement(col));
			} catch (Exception ex)
			{
			  ex.printStackTrace();
			}

			editTd.adjustWellPath(well, row, newval, wellModel);
		}
	}
}
