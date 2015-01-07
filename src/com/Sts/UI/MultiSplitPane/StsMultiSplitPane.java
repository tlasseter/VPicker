package com.Sts.UI.MultiSplitPane;

import com.Sts.UI.Beans.*;
import com.Sts.UI.MultiSplitPane.StsMultiSplitPaneLayout.*;
import com.Sts.Utilities.*;
import com.Sts.MVC.*;

import javax.accessibility.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p/>
  * All properties in this class are bound: when a properties value
  * is changed, all PropertyChangeListeners are fired.
  *
  * @author Hans Muller
  */
 public class StsMultiSplitPane extends StsJPanel
 {
     private AccessibleContext accessibleContext = null;
     private boolean continuousLayout = true;
     private DividerPainter dividerPainter = new DefaultDividerPainter();

     private boolean dragUnderway = false;
     private StsMultiSplitPaneLayout.Divider dragDivider = null;
     private Rectangle initialDividerBounds = null;
     // private boolean oldFloatingDividers = true;
     private int dragOffsetX = 0;
     private int dragOffsetY = 0;
     private int dragMin = -1;
     private int dragMax = -1;
     /**
      * Creates a StsMultiSplitPane with it's LayoutManager set to
      * to an empty MultiSplitLayout.
      */
     public StsMultiSplitPane()
     {
         super(new StsMultiSplitPaneLayout());
         InputHandler inputHandler = new InputHandler();
         addMouseListener(inputHandler);
         addMouseMotionListener(inputHandler);
         addKeyListener(inputHandler);
         setFocusable(true);
     }

	 public StsMultiSplitPane(StsMultiSplitPaneLayout l)
	 {
		 super(l);
		 InputHandler inputHandler = new InputHandler();
		 addMouseListener(inputHandler);
		 addMouseMotionListener(inputHandler);
		 addKeyListener(inputHandler);
		 setFocusable(true);
     }
     /**
      * A convenience method that returns the layout manager cast
      * to MutliSplitLayout.
      *
      * @return this StsMultiSplitPane's layout manager
      * @see java.awt.Container#getLayout
      */
     public final StsMultiSplitPaneLayout getMultiSplitLayout()
     {
         return (StsMultiSplitPaneLayout) getLayout();
     }

     public void addRootRowSplit(String name, Component component)
     {
         StsMultiSplitPaneLayout.Leaf leaf = getMultiSplitLayout().addRootRowSplit(name);
     }

	 public Component add(Component component)
	 {
		 addRootRowSplit(component.toString(),component);
		 return component;
	 }
     /**
      * A convenience method that sets the MultiSplitLayout dividerSize
      * property. Equivalent to
      * <code>getMultiSplitLayout().setDividerSize(newDividerSize)</code>.
      *
      * @param dividerSize the value of the dividerSize property
      * @see #getMultiSplitLayout
      */
     public final void setDividerSize(int dividerSize)
     {
         getMultiSplitLayout().setDividerSize(dividerSize);
     }

     /**
      * Sets the value of the <code>continuousLayout</code> property.
      * If true, then the layout is revalidated continuously while
      * a divider is being moved.  The default value of this property
      * is true.
      *
      * @param continuousLayout value of the continuousLayout property
      * @see #isContinuousLayout
      */
     public void setContinuousLayout(boolean continuousLayout)
     {
         boolean oldContinuousLayout = continuousLayout;
         this.continuousLayout = continuousLayout;
         firePropertyChange("continuousLayout", oldContinuousLayout, continuousLayout);
     }

     /**
      * Returns true if dragging a divider only updates
      * the layout when the drag gesture ends (typically, when the
      * mouse button is released).
      *
      * @return the value of the <code>continuousLayout</code> property
      * @see #setContinuousLayout
      */
     public boolean isContinuousLayout()
     {
         return continuousLayout;
     }

     /**
      * Returns the Divider that's currently being moved, typically
      * because the user is dragging it, or null.
      *
      * @return the Divider that's being moved or null.
      */
     public Divider activeDivider()
     {
         return dragDivider;
     }

     /**
      * Draws a single Divider.  Typically used to specialize the
      * way the active Divider is painted.
      */
     public static abstract class DividerPainter
     {
         /**
          * Paint a single Divider.
          *
          * @param g       the Graphics object to paint with
          * @param divider the Divider to paint
          */
         public abstract void paint(Graphics g, Divider divider);
     }

     private class DefaultDividerPainter extends DividerPainter
     {
         public void paint(Graphics g, Divider divider)
         {
             //if ((divider == activeDivider()))// && !isContinuousLayout())
             {
                 Graphics2D g2d = (Graphics2D) g;
                 g2d.setColor(Color.lightGray);
                 g2d.fill(divider.getBounds());
             }
         }
     }

     /**
      * The DividerPainter that's used to paint Dividers on this StsMultiSplitPane.
      * This property may be null.
      *
      * @return the value of the dividerPainter Property
      * @see #setDividerPainter
      */
     public DividerPainter getDividerPainter()
     {
         return dividerPainter;
     }

     /**
      * Sets the DividerPainter that's used to paint Dividers on this
      * StsMultiSplitPane.  The default DividerPainter only draws
      * the activeDivider (if there is one) and then, only if
      * continuousLayout is false.  The value of this property is
      * used by the paintChildren method: Dividers are painted after
      * the StsMultiSplitPane's children have been rendered so that
      * the activeDivider can appear "on top of" the children.
      *
      * @param dividerPainter the value of the dividerPainter property, can be null
      * @see #paintChildren
      * @see #activeDivider
      */
     public void setDividerPainter(DividerPainter dividerPainter)
     {
         this.dividerPainter = dividerPainter;
     }

     /**
      * Uses the DividerPainter (if any) to paint each Divider that
      * overlaps the clip Rectangle.  This is done after the call to
      * <code>super.paintChildren()</code> so that Dividers can be
      * rendered "on top of" the children.
      * <p/>
      * {@inheritDoc}
      */
     protected void paintChildren(Graphics g)
     {
      //   setBackground(StsModel.getCurrentModel().project.getBackgroundColor().getColor());
         super.paintChildren(g);
         DividerPainter dp = getDividerPainter();
         Rectangle clipR = g.getClipBounds();
         if ((dp != null) && (clipR != null))
         {
             Graphics dpg = g.create();
             try
             {
                 StsMultiSplitPaneLayout msl = getMultiSplitLayout();
                 for (Divider divider : msl.dividersThatOverlap(clipR))
                 {
                     dp.paint(dpg, divider);
                 }
             }
             finally
             {
                 dpg.dispose();
             }
         }
     }

     private void startDrag(int mx, int my)
     {
         requestFocusInWindow();
         StsMultiSplitPaneLayout msl = getMultiSplitLayout();
         StsMultiSplitPaneLayout.Divider divider = msl.dividerAt(mx, my);
         if (divider != null)
         {
             StsMultiSplitPaneLayout.Node prevNode = divider.previousSibling();
             StsMultiSplitPaneLayout.Node nextNode = divider.nextSibling();
             if ((prevNode == null) || (nextNode == null))
             {
                 dragUnderway = false;
             }
             else
             {
                 initialDividerBounds = divider.getBounds();
                 dragOffsetX = mx - initialDividerBounds.x;
                 dragOffsetY = my - initialDividerBounds.y;
                 dragDivider = divider;
                 Rectangle prevNodeBounds = prevNode.getBounds();
                 Rectangle nextNodeBounds = nextNode.getBounds();
                 if (dragDivider.isVertical())
                 {
                     dragMin = prevNodeBounds.x;
                     dragMax = nextNodeBounds.x + nextNodeBounds.width;
                     dragMax -= dragDivider.getBounds().width;
                 }
                 else
                 {
                     dragMin = prevNodeBounds.y;
                     dragMax = nextNodeBounds.y + nextNodeBounds.height;
                     dragMax -= dragDivider.getBounds().height;
                 }
                 // oldFloatingDividers = getMultiSplitLayout().getFloatingDividers();
                 // StsException.systemDebug(this, "startDrag", "oldFloatingDividers: set to " + oldFloatingDividers);
                 // getMultiSplitLayout().setFloatingDividers(false);
                 getMultiSplitLayout().setDraggingDivider(true);
                 dragUnderway = true;
             }
         }
         else
         {
             dragUnderway = false;
         }
     }

     private void repaintDragLimits()
     {
         Rectangle damageR = dragDivider.getBounds();
         if (dragDivider.isVertical())
         {
             damageR.x = dragMin;
             damageR.width = dragMax - dragMin;
         }
         else
         {
             damageR.y = dragMin;
             damageR.height = dragMax - dragMin;
         }
         repaint(damageR);
     }

     private void updateDrag(int mx, int my)
     {
         if (!dragUnderway)
         {
             return;
         }
         Rectangle oldBounds = dragDivider.getBounds();
         Rectangle bounds = new Rectangle(oldBounds);
         if (dragDivider.isVertical())
         {
             bounds.x = mx - dragOffsetX;
             bounds.x = Math.max(bounds.x, dragMin);
             bounds.x = Math.min(bounds.x, dragMax);
         }
         else
         {
             bounds.y = my - dragOffsetY;
             bounds.y = Math.max(bounds.y, dragMin);
             bounds.y = Math.min(bounds.y, dragMax);
         }
         dragDivider.setBounds(bounds);
         // if(StsMultiSplitPaneLayout.debug) StsException.systemDebug(this, "updateDrag", " dragDivider.bounds.x : " + bounds.x);
         if (isContinuousLayout())
         {
             revalidate();
             repaintDragLimits();
         }
         else
         {
             repaint(oldBounds.union(bounds));
         }
     }

     private void clearDragState()
     {
         dragDivider = null;
         initialDividerBounds = null;
         // oldFloatingDividers = true;
         // StsException.systemDebug(this, "clearDragState", "oldFloatingDividers: set to true");
         dragOffsetX = dragOffsetY = 0;
         dragMin = dragMax = -1;
         dragUnderway = false;
         getMultiSplitLayout().setDraggingDivider(false);
     }

     private void finishDrag(int x, int y)
     {
         if (dragUnderway)
         {
             clearDragState();
             if (!isContinuousLayout())
             {
                 revalidate();
                 repaint();
             }
         }
     }

     private void cancelDrag()
     {
         if (dragUnderway)
         {
             clearDragState();
             dragDivider.setBounds(initialDividerBounds);
             // getMultiSplitLayout().setFloatingDividers(oldFloatingDividers);
             getMultiSplitLayout().setDraggingDivider(false);
             setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
             repaint();
             revalidate();
         }
     }

     private void updateCursor(int x, int y, boolean show)
     {
         if (dragUnderway)
         {
             return;
         }
         int cursorID = Cursor.DEFAULT_CURSOR;
         if (show)
         {
             StsMultiSplitPaneLayout.Divider divider = getMultiSplitLayout().dividerAt(x, y);
             if (divider != null)
             {
                 cursorID = (divider.isVertical()) ?
                     Cursor.E_RESIZE_CURSOR :
                     Cursor.N_RESIZE_CURSOR;
             }
         }
         setCursor(Cursor.getPredefinedCursor(cursorID));
     }


     private class InputHandler extends MouseInputAdapter implements KeyListener
     {

         public void mouseEntered(MouseEvent e)
         {
             updateCursor(e.getX(), e.getY(), true);
         }

         public void mouseMoved(MouseEvent e)
         {
             updateCursor(e.getX(), e.getY(), true);
         }

         public void mouseExited(MouseEvent e)
         {
             updateCursor(e.getX(), e.getY(), false);
         }

         public void mousePressed(MouseEvent e)
         {
             startDrag(e.getX(), e.getY());
         }

         public void mouseReleased(MouseEvent e)
         {
             finishDrag(e.getX(), e.getY());
         }

         public void mouseDragged(MouseEvent e)
         {
             updateDrag(e.getX(), e.getY());
         }

         public void keyPressed(KeyEvent e)
         {
             if (e.getKeyCode() == KeyEvent.VK_ESCAPE)
             {
                 cancelDrag();
             }
         }

         public void keyReleased(KeyEvent e) { }

         public void keyTyped(KeyEvent e) { }
     }

     public AccessibleContext getAccessibleContext()
     {
         if (accessibleContext == null)
         {
             accessibleContext = new AccessibleMultiSplitPane();
         }
         return accessibleContext;
     }

     protected class AccessibleMultiSplitPane extends AccessibleJPanel
     {
         public AccessibleRole getAccessibleRole()
         {
             return AccessibleRole.SPLIT_PANE;
         }
     }

     public void removeLeaf(String leafName)
     {
         StsMultiSplitPaneLayout.Leaf leaf = getLeaf(leafName);
         if (leaf == null)
         {
             StsException.systemError(this, "removeLeaf", "Can't find leaf to remove: " + leafName);
             return;
         }
         if (leaf.isSplit())
         {
             StsException.systemError(this, "removeLeaf", "Can't remove leaf " + leafName + ". already split.");
             return;
         }
         Component leafComponent = leaf.getComponent();

         remove(leafComponent);
         StsMultiSplitPaneLayout layout = getMultiSplitLayout();
         layout.removeLayoutLeaf(leafName);
     }

     public StsMultiSplitPaneLayout.Leaf getLeaf(String leafName)
     {
         StsMultiSplitPaneLayout layout = getMultiSplitLayout();
         return layout.getLeaf(leafName);
     }

     public void removeLeafComponent(String leafName)
     {
         StsMultiSplitPaneLayout layout = getMultiSplitLayout();
         StsMultiSplitPaneLayout.Leaf leaf = layout.getLeaf(leafName);
         Component leafComponent = leaf.getComponent();
         remove(leafComponent);
         leaf.removeComponent();
     }

     public void replaceLeafComponent(String leafName, Component component)
     {
         StsMultiSplitPaneLayout layout = getMultiSplitLayout();
         StsMultiSplitPaneLayout.Leaf leaf = layout.getLeaf(leafName);
         Component leafComponent = leaf.getComponent();
         if(leafComponent == component) return;
         remove(leafComponent);
         leaf.setComponent(component);
         add(component, leafName);
     }

     public void addLeafComponent(String leafName, Component component)
     {
         StsMultiSplitPaneLayout layout = getMultiSplitLayout();
         StsMultiSplitPaneLayout.Leaf leaf = layout.getLeaf(leafName);
         leaf.setComponent(component);
         add(component, leafName);
     }

     public void addRootLeaf(String leafName)
     {
         getMultiSplitLayout().addRootLeaf(leafName);
     }

     public void addRootLeafComponent(String leafName, Component component)
     {
         getMultiSplitLayout().addRootLeaf(leafName);
         add(component, leafName);
     }

     static public void main(String[] args)
     {
         final StsMultiSplitPane multiSplitPane = new StsMultiSplitPane();
         final StsMultiSplitPaneLayout layout = multiSplitPane.getMultiSplitLayout();
         final Dimension windowSize = new Dimension(200, 200);
         final Dimension leftPanelSize = new Dimension(200, 200);
         final Dimension rightPanelSize = new Dimension(200, 200);
         multiSplitPane.setSize(windowSize);
		 //multiSplitPane.setColor(Color.WHITE);
         Runnable runnable = new Runnable()
         {
             public void run()
             {
                 layout.addRootLeaf("panel");
                 //               layout.addRootRowSplit("left");
                 StsJPanel panelLeft = new StsJPanel();
                 JButton leftButton = new JButton("left");
                 panelLeft.add(leftButton);
                 leftButton.setPreferredSize(leftPanelSize);
                 multiSplitPane.add(panelLeft, "panel");
				 layout.addRootRowSplit("left2");
				 StsJPanel panelLeft2 = new StsJPanel();
				 JButton left2Button = new JButton("left2");
				 panelLeft2.add(left2Button);
				 left2Button.setPreferredSize(leftPanelSize);
				 multiSplitPane.add(panelLeft2, "left2");


             }
         };
         StsToolkit.runWaitOnEventThread(runnable);

         final JDialog dialog = StsToolkit.createDialog(multiSplitPane, false);

         try
         {
             Thread.sleep(3000);
         }
         catch (Exception e)
         {

         }

         runnable = new Runnable()
         {
             public void run()
             {
				 System.out.println("add right");
                 layout.addRootRowSplit("right");
                 StsJPanel panelRight = new StsJPanel();
                 JButton rightButton = new JButton("right2");
                 panelRight.add(rightButton);
				 rightButton.setPreferredSize(rightPanelSize);
				 JButton rightButton2 = new JButton("right3");
                 panelRight.add(rightButton2);
                 rightButton2.setPreferredSize(rightPanelSize);
                 multiSplitPane.add(panelRight, "right");
                 dialog.pack();
                 dialog.repaint();
				 //layout.removeLayoutLeaf("right");
				 multiSplitPane.remove(panelRight);
             }
         };
         StsToolkit.runWaitOnEventThread(runnable);

     }
     /*
 static public void main(String[] args)
     {
         final StsMultiSplitPane multiSplitPane = new StsMultiSplitPane();
         final StsMultiSplitPaneLayout layout = multiSplitPane.getMultiSplitLayout();

         Runnable runnable = new Runnable()
         {
             public void run()
             {
                 layout.addRootRowSplit(new String[] { "left", "center"} );
                 multiSplitPane.add(new JButton("Left Component"), "left");
                 multiSplitPane.add(new JButton("Center Component"), "center");
             }
         };
         StsToolkit.runWaitOnEventThread(runnable);

         final JDialog dialog = StsToolkit.createDialog(multiSplitPane, false);
         try
         {
             Thread.sleep(3000);
         }
         catch(Exception e)
         {

         }

         runnable = new Runnable()
         {
             public void run()
             {
                 layout.addRootChild("right");
                 multiSplitPane.add(new JButton("Right Component"), "right");
                 dialog.pack();
 //                dialog.invalidate();
 //                dialog.validate();
                 dialog.repaint();
             }
         };
         StsToolkit.runWaitOnEventThread(runnable);
     }
     */

/* static public void main(String[] args)
     {
         final StsMultiSplitPane multiSplitPane = new StsMultiSplitPane();
         final StsMultiSplitPaneLayout layout = multiSplitPane.getMultiSplitLayout();
         layout.addRootRowSplit(new String[] { "left", "center", "right"} );
         layout.addColSplit("center", new String[] { "top", "bottom"});

         Runnable runnable = new Runnable()
         {
             public void run()
             {
                 multiSplitPane.add(new JButton("Left Component"), "left");
                 multiSplitPane.add(new JButton("Right Component"), "right");
                 multiSplitPane.add(new JButton("Center Component"), "center");
                 multiSplitPane.add(new JButton("Top Component"), "top");
                 multiSplitPane.add(new JButton("Bottom Component"), "bottom");
             }
         };
         StsToolkit.runWaitOnEventThread(runnable);

         final JDialog dialog = StsToolkit.createDialog(multiSplitPane, false);
         try
         {
             Thread.sleep(3000);
         }
         catch(Exception e)
         {

         }
         runnable = new Runnable()
         {
             public void run()
             {
 //                layout.removeLayoutLeaf("center");
                 multiSplitPane.removeLeaf("top");
                 dialog.pack();
                 dialog.invalidate();
                 dialog.validate();
                 dialog.repaint();
             }
         };
         StsToolkit.runWaitOnEventThread(runnable);
     }
 */
 }
