package com.Sts.UI;

/**
 * Created by IntelliJ IDEA.
 * User: S2S Systems
 * Date: Sep 2, 2011
 * Time: 11:34:48 AM
 * To change this template use File | Settings | File Templates.
 */
import com.Sts.DBTypes.StsMainTimeObject;
import com.Sts.DBTypes.StsObject;
import com.Sts.MVC.StsModel;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Enumeration;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.table.*;

public class StsTimeTable extends JTable
{
  private final int ROWDISTANCE = 0;    // the distance between the rows
  private final int DATEROWHEIGHT = 16;   // height of the date row
  private ImageTableModel itm;      // table model
  private DefaultTableColumnModel colMdl = new DefaultTableColumnModel();
  private int maxRowHeight = 16;      // maximum row height
  private boolean autoSizeRowHeight = false;    // resize row height automatically
  private StsModel model = null;
  /**
   * Generates a new TimeLine object
   * @param rows  - number of rows of the default grid
   * @param columns - number of columns of the default grid
   * @param autoSizeRowHeight  - auto adjust row height to the max row height
   */
  public StsTimeTable(StsModel mdl, int rows, int columns, boolean autoSizeRowHeight){
    this.autoSizeRowHeight = autoSizeRowHeight;
    this.model= mdl;
    setTimeLineLookAndFeel();
    itm = initModel();
    itm.initGrid(rows, columns);
    setModel(itm);  // set the table model
  }
  /**
   * Adds a new column to the table
   * @param date - The column name as java.util.Date
   * @param columndata - The row values for this column.
   * Object[0] contains the cell text, object[1] contains the ImageIcon
   */
  public void addColumn(java.util.Date date, boolean includeDate, boolean includeTime, ArrayList<Object[]> columndata){
    itm.addColumn(convertDate(date, includeDate, includeTime), columndata);
    // adjust row height
    //if(autoSizeRowHeight){
    //  calcMaxRowHeight(columndata);
    //}
  }

  /**
   * Adds a new column to the table
   * @param columnName - The column name
   * @param columndata - The row values for this column.
   * Object[0] contains the cell text, object[1] contains the ImageIcon
   */
  public void addColumn(String columnName, ArrayList<Object[]> columndata){
    itm.addColumn(columnName, columndata);
    // adjust row height
    //if(autoSizeRowHeight){
    //  calcMaxRowHeight(columndata);
    //}
  }

    public void setColumnWidths()
    {
        Enumeration cols = getColumnModel().getColumns();
        while(cols.hasMoreElements())
        {
            TableColumn col = (TableColumn)cols.nextElement();
            col.setMinWidth(5);
            col.setMaxWidth(5);
            col.setWidth(5);
            col.setPreferredWidth(5);
        }
    }
  /**
   * Utility function to set the scroll pane
   */
  public int getRowCount(){
    return itm.getRowCount();
  }

  /**
   * Utility function to set the scroll pane
   */
  public int getColumnCount(){
    return itm.getColumnCount();
  }

  /**
   * @return ImageTableModel - a new ImageTableModel
   */
  private ImageTableModel initModel(){
    return new ImageTableModel();
  }
  /**
   * Set some JTable properties to make it
   * look more like a timeline
   */
  private void setTimeLineLookAndFeel(){
    this.getTableHeader().setReorderingAllowed(false);
    this.setCellSelectionEnabled(true);
    this.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    this.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    //this.setIntercellSpacing(new Dimension(0,0));
    this.setShowGrid(true);
    this.setGridColor(new Color(200,200,200));
    this.setShowHorizontalLines(true);
    this.setTableHeader(null);
    this.setDefaultRenderer(Object.class, new MyCellRenderer(this)); 
    this.addMouseMotionListener(new MyMouseMotionAdapter(this));
  }
  /**
   * Converts a java.util.Date object to
   * a readable format
   * @param date - a java.util.Date object
   * @return String - Date as a String
   */
  public String convertDate(java.util.Date date, boolean includeDate, boolean includeTime){
    SimpleDateFormat df = null;
    if(includeDate && includeTime)
         df = new SimpleDateFormat("dd/MM/yy hh:mm:ss");
    else if(includeDate && !includeTime)
        df = new SimpleDateFormat("dd/MM/yy");
    else if(!includeDate && includeTime)
        df = new SimpleDateFormat("hh:mm:ss");
    else
        df = new SimpleDateFormat("dd/MM/yy");
    return df.format(date);
  }
  /**
   * Calculates the maximum row height
   * @param list - the column list
   */
  private void calcMaxRowHeight(ArrayList<Object[]> list){
    int max = 0;
    for(int i=0;i<list.size();i++){
      int labelheight = 0;
      int imageheight = 0;
      if(list.get(i) != null){
        if(list.get(i)[0] != null){
          if(list.get(i)[0] instanceof String){
            labelheight = new JLabel((String)list.get(i)[0]).getPreferredSize().height;
          }
        }
        if(list.get(i)[1] != null){
          if(list.get(i)[1] instanceof ImageIcon){
            imageheight = ((ImageIcon)list.get(i)[1]).getIconHeight();
          }
        }
      }
      int height = labelheight + imageheight + ROWDISTANCE;
      if(maxRowHeight < height){
        maxRowHeight = height;
      }
    }
  }
  @Override
  public void paintComponent(Graphics g){
    super.paintComponent(g);
    g.setColor( Color.BLACK );
    int y = getRowHeight(0)*(itm.getRowCount()-1) - 1;
    g.drawLine(0, y, getSize().width, y);
  }
  //@Override
  /*
  public Component prepareRenderer(TableCellRenderer renderer, int row, int column){
    Component c = super.prepareRenderer(renderer, row, column);
    if(row == itm.getRowCount()-1){
      c.setBackground(java.awt.SystemColor.control);
    }
    else if( !this.isCellSelected(row, column)){
      c.setBackground(column % 2 != 0 ? new Color(241, 245, 250) : null);
    }
    return c;
  }
    */
  @Override
  /**
   * makes the date column look like a table header
   */
  public void changeSelection(int row, int column, boolean toggle, boolean extend){
    Object[] o = (Object[])getValueAt(row, column);
    if (o != null && row != itm.getRowCount()-1){
      super.changeSelection(row, column, toggle, extend);
    }
    if(o.length != 3)
        return;
    StsMainTimeObject obj = (StsMainTimeObject)model.getObjectWithName(model.getStsClass("com.Sts.DBTypes.Sts" + (String)o[2]), (String)o[0]);
    model.getProject().setProjectTime(obj.getBornDateLong(), true);
  }

  /**
   *   MouseMotionAdapter for this table
   *   Used to set the mouse cursor only
   */
  class MyMouseMotionAdapter extends MouseMotionAdapter{
    private StsTimeTable theTimeLine;
    public MyMouseMotionAdapter(StsTimeTable table){
      this.theTimeLine = table;
    }
    public void mouseMoved(MouseEvent me) {
      int row = rowAtPoint( me.getPoint() );
      if(row < theTimeLine.itm.getRowCount()-1){
        theTimeLine.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
      }
      else{
        theTimeLine.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
    }
  }

  /**
   *  Cell renderer for this table
   */
  class MyCellRenderer extends DefaultTableCellRenderer{
    private StsTimeTable theTimeLine;
    public MyCellRenderer(StsTimeTable timeline){
      this.theTimeLine = timeline;
      setHorizontalAlignment( CENTER );
      setVerticalTextPosition(JLabel.BOTTOM);
      setHorizontalTextPosition(JLabel.CENTER);
      setFont(new java.awt.Font("Dialog", 0, 8));
    }

    public void setValue(Object value)
    {
      Object[] a = (Object[])value;
      if(value == null)
      {      // empty cell
        setIcon(null);
        super.setValue(value);
      }
      else if (a[1] instanceof Icon) {    // normal cell with text and image
        setIcon((Icon) a[1]);
        //setText((String)a[0]);
        //this.setAlignmentY(JLabel.BOTTOM_ALIGNMENT);
        //theTimeLine.setRowHeight(maxRowHeight);
        //theTimeLine.setRowHeight(itm.getRowCount()-1,DATEROWHEIGHT);
        setText(null);
        theTimeLine.setToolTipText((String)a[0]);
        //this.setBorder(new RoundBorder(Color.red,Color.green));
      }
      else if(a[1] == null && a[0] instanceof String){
        setText((String)a[0]);
        setIcon(null);
      }
      else {
        setIcon(null);
        super.setValue(value);
      }
    }
  }

  /**
   * The table model for this timeline   *
   */
  class ImageTableModel extends AbstractTableModel{
    private ArrayList<String> columnnames; // holds the column names
    private ArrayList<ArrayList<Object[]>> data; // holds the table data
    private int maxRowCount;
    private int columnCursor; // points on the current column
    public ImageTableModel(){
      columnnames = new ArrayList<String>();
      data = new ArrayList<ArrayList<Object[]>>();
      maxRowCount = 0;
      columnCursor = 0;
      autoResizeMode = AUTO_RESIZE_OFF;
    }
    public Object getValueAt(int row, int column){
      if (data.get(column).size()-1<row){
        return null;
      }
      else{
        return data.get(column).get(row);
      }
    }
    public int getRowCount(){
      return maxRowCount;
    }
    public int getColumnCount(){
      return columnnames.size();
    }
    public String getColumnName( int columnIndex ){
      return columnnames.get(columnIndex);
    }
    /**
     * Adds a new column to the table
     * @param columnName - The column name
     * @param columndata - The row values for this column.
     */
    public void addColumn(String columnName, ArrayList<Object[]> columndata){
      if(columnCursor >= columnnames.size()){
        columnnames.add(columnName);
        data.add(rotateFillList(columnName,columndata));
      }
      else{
        columnnames.set(columnCursor, columnName);
        data.set(columnCursor, rotateFillList(columnName,columndata));
      }
      SwingUtilities.invokeLater (new Runnable(){    // fixes a nasty java vector bug
        public void run ()  {
          fireTableStructureChanged();
          fireTableDataChanged();
        }
      });
      columnCursor++;
    }
    public void initGrid(int rows, int columns){
      for(int i=0;i<columns;i++){
        ArrayList<Object[]> newdata = new ArrayList<Object[]>();
        for(int j=0;j<rows;j++){
          newdata.add(null);
        }
        columnnames.add(String.valueOf(i));
        data.add(newdata);
        maxRowCount = rows;
      }
      SwingUtilities.invokeLater (new Runnable(){    // fixes a nasty java vector bug
        public void run ()  {
          fireTableStructureChanged();
          fireTableDataChanged();
        }
      });
    }
    /**
     * Rotates the list. If list.size() is smaller than
     * maxRowCount the list if filled with null values
     * This generates the bottom up effect
     * @param columnName - The column name
     * @param list
     * @return list
     */
    private ArrayList<Object[]> rotateFillList(String columnName, ArrayList<Object[]> list){
      if(columnName != null)
          list.add(0,new Object[]{columnName,null});  // set column name to be on the bottom
      if(maxRowCount < list.size()){
        // adjust all rows to the new maxRowCount
        maxRowCount = list.size();
        for(int i=0;i<data.size();i++){
          int diff = maxRowCount - data.get(i).size();
          for(int j=0;j<diff;j++){
            data.get(i).add(0,null);
          }
        }
      }
      else {  // fill with null values
        int diff = maxRowCount - list.size();
        for(int i=0;i<diff;i++){
          list.add(null);
        }
      }
      ArrayList<Object[]> rotatedList = new ArrayList<Object[]>();
      for(int i= list.size()-1;i>=0;i--){    // rotate list
        rotatedList.add(list.get(i));
      }
      return rotatedList;
    }
  }

}
