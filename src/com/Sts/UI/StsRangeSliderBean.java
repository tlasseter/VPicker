//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.Utilities.StsMath;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
 
public class StsRangeSliderBean extends JComponent implements MouseListener, MouseMotionListener, KeyListener {

         protected class StsRangeSliderChangeListener implements ChangeListener {
                 public void stateChanged(ChangeEvent e) {
                         fireChangeEvent();
                 }
         }
         protected final static int INT = 0;
         protected final static int FLOAT = 1;
         
         protected final static int ARROW_HEIGHT = 4;
         protected final static int ARROW_SZ = 16;
         protected final static int ARROW_WIDTH = 8;
         private static final BasicStroke arrowStroke = new BasicStroke(0);
         public static final int HORIZONTAL = 1;
         public static final int LEFTRIGHT_TOPBOTTOM = 0;
         private static final int PICK_LEFT_OR_TOP = 1;
         private static final int PICK_NONE = 0;
         private static final int PICK_RIGHT_OR_BOTTOM = 3;
         private static final int PICK_THUMB = 2;
         public static final int PREFERRED_BREADTH = 16;
         public static final int PREFERRED_LENGTH = 300;
         public static final int RIGHTLEFT_BOTTOMTOP = 1;
         public static final int VERTICAL = 0;
         private Color arrowColor = new Color(77, 97, 133);
         protected ChangeEvent changeEvent = null;
         protected ChangeListener changeListener;
         protected int direction;
         protected int numType = INT;
         protected float numScale = 10.0f;
         protected boolean empty;
         protected int increment = 1;
         protected ArrayList listeners = new ArrayList();
         protected int minExtent = 0; // min extent, in pixels
         protected BoundedRangeModel model;
         int mouse;
         protected int orientation;
         int pick;
         int pickOffsetHigh;
         int pickOffsetLow;
         protected Color thumbColor = new Color(150, 180, 220);
         private boolean use3DHighlights = false;

         public StsRangeSliderBean(BoundedRangeModel model, int orientation, int direction) 
         {
                 super.setFocusable(true);
                 this.model = model;
                 this.orientation = orientation;
                 this.direction = direction;
                 setForeground(Color.LIGHT_GRAY);
 
                 this.changeListener = createListener();
                 model.addChangeListener(changeListener);
 
                 addMouseListener(this);
                 addMouseMotionListener(this);
                 addKeyListener(this);
         }

         // Since the range model is integer, this constructor will determine the best scalar to mulitple the float min/max by to
         // maintain resolution.
         public StsRangeSliderBean(float minimum, float maximum, float lowValue, boolean autoCompute, int orientation)
         {
             float raise = 10.0f;
             if(autoCompute)
             {
                int sigDigits = (int) StsMath.max(new float[] {Math.abs(maximum), Math.abs(minimum)});
                String remainStg = String.valueOf(sigDigits);
                int power = 9-remainStg.length();
                if(sigDigits == 0)
                    power = 9;
                raise = (float)Math.pow((double)10,(double)power);
             }
             lowValue = lowValue * raise;
             minimum = minimum * raise;
             maximum = maximum * raise;

             super.setFocusable(true);
             this.model = new DefaultBoundedRangeModel((int)lowValue, (int)maximum-(int)minimum, (int)minimum, (int)maximum);
             this.orientation = orientation;
             this.direction = LEFTRIGHT_TOPBOTTOM;
             setForeground(Color.LIGHT_GRAY);

             this.changeListener = createListener();
             model.addChangeListener(changeListener);

             addMouseListener(this);
             addMouseMotionListener(this);
             addKeyListener(this);
        	 numType = FLOAT;
             numScale = raise;
         }

         public StsRangeSliderBean(float minimum, float maximum, float lowValue, int orientation)
         {
             this(new DefaultBoundedRangeModel((int)(lowValue*10.0f), (int)(maximum*10.0f)-(int)(minimum*10.0f), (int)(minimum*10.0f), (int)(maximum*10.0f)), orientation,
                             LEFTRIGHT_TOPBOTTOM);
        	 numType = FLOAT;
             numScale = 10.0f;
         } 
         
         public StsRangeSliderBean(int minimum, int maximum, int lowValue, int orientation) {
                 this(new DefaultBoundedRangeModel(lowValue, maximum-minimum, minimum, maximum), orientation,
                                 LEFTRIGHT_TOPBOTTOM);
         }
 
         public StsRangeSliderBean(int minimum, int maximum, int lowValue, int orientation, int direction) {
                 this(new DefaultBoundedRangeModel(lowValue, maximum-minimum, minimum, maximum), orientation, direction);
         }
 
         public void addChangeListener(ChangeListener cl) {
                 if (!listeners.contains(cl))
                         listeners.add(cl);
         }
 
         protected ChangeListener createListener() {
                 return new StsRangeSliderChangeListener();
         }
 
         protected void customPaint(Graphics2D g, int width, int height) {
                 // does nothing in this class
                 // subclasses can override to perform custom painting
         }

         public float getScalar() { return numScale; }
         protected void fireChangeEvent() {
                 repaint();
                 if (changeEvent == null) {
                         changeEvent = new ChangeEvent(this);
                 }
                 Iterator iter = listeners.iterator();
                 while (iter.hasNext()) {
                         ((ChangeListener) iter.next()).stateChanged(changeEvent);
                 }
         }
         
         public float getFloatHighValue() 
         {
             return (float)(model.getValue()/numScale) + (float)(model.getExtent()/numScale);
         }
 
         public float getFloatLowValue() 
         {
                 return (float)(model.getValue()/numScale);
         } 
         
         public int getHighValue() 
         {
             return model.getValue() + model.getExtent();
         }
 
         public int getLowValue() {
                 return model.getValue();
         }
         
         public float getFloatMaximum() {
             return (float)(model.getMaximum()/numScale);
         }

         public float getFloatMinimum() {
             return (float)(model.getMinimum()/numScale);
         } 

         public void setFloatMinimum(float min)
         {
             int scaledMin = (int)(min * numScale);
             model.setMinimum(scaledMin);
             model.setExtent(model.getMaximum()-model.getMinimum());
             model.setValue(model.getMinimum());
         }

         public int getMaximum() {
                 return model.getMaximum();
         }

         public void setFloatMaximum(float max)
         {
             int scaledMax = (int)(max * numScale);
             model.setMaximum(scaledMax);
             model.setExtent(model.getMaximum()-model.getMinimum());
             model.setValue(model.getMinimum());
         }

         public int getMinimum() {
                 return model.getMinimum();
         }
 
         public BoundedRangeModel getModel() {
                 return model;
         }

         public Dimension getPreferredSize() {
                 if (orientation == VERTICAL) {
                         return new Dimension(PREFERRED_BREADTH, PREFERRED_LENGTH);
                 } else {
                         return new Dimension(PREFERRED_LENGTH, PREFERRED_BREADTH);
                 }
         }
 
         public Color getThumbColor() {
                 return thumbColor;
         }
 
         private void grow(int increment) {
                 model.setRangeProperties(model.getValue() - increment, model.getExtent() + 2 * increment, model
                                 .getMinimum(), model.getMaximum(), false);
         }
 
         // ------------------------------------------------------------------------
         // Rendering
 
         public void keyPressed(KeyEvent e) {
                 int kc = e.getKeyCode();
                 boolean v = (orientation == VERTICAL);
                 boolean d = (kc == KeyEvent.VK_DOWN);
                 boolean u = (kc == KeyEvent.VK_UP);
                 boolean l = (kc == KeyEvent.VK_LEFT);
                 boolean r = (kc == KeyEvent.VK_RIGHT);
 
                 int minimum = getMinimum();
                 int maximum = getMaximum();
                 int lowValue = getLowValue();
                 int highValue = getHighValue();
 
                 if (v && r || !v && u) {
                         if (lowValue - increment >= minimum && highValue + increment <= maximum) {
                                 grow(increment);
                         }
                 } else if (v && l || !v && d) {
                         if (highValue - lowValue >= 2 * increment) {
                                 grow(-1 * increment);
                         }
                 } else if (v && d || !v && l) {
                         if (lowValue - increment >= minimum) {
                                 offset(-increment);
                         }
                 } else if (v && u || !v && r) {
                         if (highValue + increment <= maximum) {
                                 offset(increment);
                         }
                 }
         }
 
         public void keyReleased(KeyEvent e) {
         }
 
         public void keyTyped(KeyEvent e) {
         }
 
         public void mouseClicked(MouseEvent e) {
         }
 
         public void mouseDragged(MouseEvent e) {
                 requestFocus();
                 int value = (orientation == VERTICAL) ? e.getY() : e.getX();
 
                 int minimum = getMinimum();
                 int maximum = getMaximum();
                 int lowValue = getLowValue();
                 int highValue = getHighValue();
 
                 switch (pick) {
                 case PICK_LEFT_OR_TOP:
                         int low = toLocal(value - pickOffsetLow);
 
                         if (low < minimum) {
                                 low = minimum;
                         }
                         if (low > maximum) {
                                 low = maximum;
                         }
                         if (low > highValue - minExtent) {
                                 low = highValue - minExtent;
                         }
                         setLowValue(low);
                         break;
 
                 case PICK_RIGHT_OR_BOTTOM:
                         int high = toLocal(value - pickOffsetHigh);
 
                         if (high < minimum) {
                                 high = minimum;
                         }
                         if (high > maximum) {
                                 high = maximum;
                         }
                         if (high < lowValue + minExtent) {
                                 high = lowValue + minExtent;
                         }
                         setHighValue(high);
                         break;
 
                 case PICK_THUMB:
                         int dxOrDy = toLocal(value - pickOffsetLow) - lowValue;
                         if ((dxOrDy < 0) && ((lowValue + dxOrDy) < minimum)) {
                                 dxOrDy = minimum - lowValue;
                         }
                         if ((dxOrDy > 0) && ((highValue + dxOrDy) > maximum)) {
                                 dxOrDy = maximum - highValue;
                         }
                         if (dxOrDy != 0) {
                                 offset(dxOrDy);
                         }
                         break;
                 }
         }
 
         public void mouseEntered(MouseEvent e) {
         }
 
         public void mouseExited(MouseEvent e) {
         }
 
         // ------------------------------------------------------------------------
         // Event Handling
 
         public void mouseMoved(MouseEvent e) {
                 if (orientation == VERTICAL) {
                         switch (pickHandle(e.getY())) {
                         case PICK_LEFT_OR_TOP:
                                 setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                                 break;
                         case PICK_RIGHT_OR_BOTTOM:
                                 setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                                 break;
                         case PICK_THUMB:
                                 setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                                 break;
                         case PICK_NONE:
                                 setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                                 break;
                         }
                 } else {
                         switch (pickHandle(e.getX())) {
                         case PICK_LEFT_OR_TOP:
                                 setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                                 break;
                         case PICK_RIGHT_OR_BOTTOM:
                                 setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                                 break;
                         case PICK_THUMB:
                                 setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                                 break;
                         case PICK_NONE:
                                 setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                                 break;
                         }
                 }
         }
 
         public void mousePressed(MouseEvent e) {
                 if (orientation == VERTICAL) {
                         pick = pickHandle(e.getY());
                         pickOffsetLow = e.getY() - toScreen(getLowValue());
                         pickOffsetHigh = e.getY() - toScreen(getHighValue());
                         mouse = e.getY();
                 } else {
                         pick = pickHandle(e.getX());
                         pickOffsetLow = e.getX() - toScreen(getLowValue());
                         pickOffsetHigh = e.getX() - toScreen(getHighValue());
                         mouse = e.getX();
                 }
                 repaint();
         }
 
         public void mouseReleased(MouseEvent e) {
                 pick = PICK_NONE;
                 repaint();
         }
 
         private void offset(int dxOrDy) {
                 model.setValue(model.getValue() + dxOrDy);
         }
 
         protected void paint3DRectLighting(Graphics2D g2, int x, int y, int width, int height) {
                 if (use3DHighlights) {
                         g2.setColor(Color.white);
                         g2.drawLine(x + 1, y + 1, x + 1, y + height - 1);
                         g2.drawLine(x + 1, y + 1, x + width - 1, y + 1);
                         g2.setColor(Color.gray);
                         g2.drawLine(x + 1, y + height - 1, x + width - 1, y + height - 1);
                         g2.drawLine(x + width - 1, y + 1, x + width - 1, y + height - 1);
                         g2.setColor(Color.darkGray);
                         g2.drawLine(x, y + height, x + width, y + height);
                         g2.drawLine(x + width, y, x + width, y + height);
                 }
         }
 
         protected void paintArrow(Graphics2D g2, double x, double y, int w, int h, boolean topDown) {
                 int intX = (int) (x + 0.5);
                 int intY = (int) (y + 0.5);
 
                 if (orientation == VERTICAL) {
                         if (w % 2 == 0) {
                                 w = w - 1;
                         }
 
                         if (topDown) {
                                 for (int i = 0; i < (w / 2 + 1); i++) {
                                         g2.drawLine(intX + i, intY + i, intX + w - i - 1, intY + i);
                                 }
                         } else {
                                 for (int i = 0; i < (w / 2 + 1); i++) {

                                         g2.drawLine(intX + w / 2 - i, intY + i, intX + w - w / 2 + i - 1, intY + i);

                                 }
                         }
                 } else {
                         if (h % 2 == 0) {
                                 h = h - 1;
                         }
 
                         // do this for vertical also (ronyeh)
                         final int x_1 = intX + 1;
                         final int h_1 = h + 1;
                         final int h_over_2_plus_1 = Math.round((float)(h / 2.0)) + 1;
                         final int h_1_5_minus_1 = Math.round((float)(1.5 * h)) - 1;
                         final int x_w_minus_1 = intX + w - 1;
                         final int x_w_minus_2 = intX + w - 2;
 
                         if (topDown) {
                                 // leftmost
                                 g2.drawLine(intX, h_over_2_plus_1, x_w_minus_2, h);
                                 g2.drawLine(intX, h_1_5_minus_1, x_w_minus_2, h_1);
 
                                 // middle
                                 g2.drawLine(x_1, h_over_2_plus_1, x_w_minus_1, h);
                                 g2.drawLine(x_1, h_1_5_minus_1, x_w_minus_1, h_1);
 
                                 // right most
                                 g2.drawLine(x_1, h_over_2_plus_1 - 1, intX + w, h);
                                 g2.drawLine(x_1, h_1_5_minus_1 + 1, intX + w, h_1);
                         } else {
                                 g2.drawLine(x_w_minus_1, h_over_2_plus_1, x_1, h);
                                 g2.drawLine(x_w_minus_1, h_1_5_minus_1, x_1, h_1);
 
                                 // middle
                                 g2.drawLine(x_w_minus_2, h_over_2_plus_1, intX, h);
                                 g2.drawLine(x_w_minus_2, h_1_5_minus_1, intX, h_1);
 
                                 // left most
                                 g2.drawLine(x_w_minus_2, h_over_2_plus_1 - 1, intX - 1, h);
                                 g2.drawLine(x_w_minus_2, h_1_5_minus_1 + 1, intX - 1, h_1);
                         }
                 }
         }
 
         public void paintComponent(Graphics g) {
                 Rectangle bounds = getBounds();
                 int width = (int) bounds.getWidth() - 1;
                 int height = (int) bounds.getHeight() - 1;
 
                 int min = toScreen(getLowValue());
                 int max = toScreen(getHighValue());
 
                 // Paint the full slider if the slider is marked as empty
                 if (empty) {
                         if (direction == LEFTRIGHT_TOPBOTTOM) {
                                 min = ARROW_SZ;
                                 max = (orientation == VERTICAL) ? height - ARROW_SZ : width - ARROW_SZ;
                         } else {
                                 min = (orientation == VERTICAL) ? height - ARROW_SZ : width - ARROW_SZ;
                                 max = ARROW_SZ;
                         }
                 }
 
                 final Graphics2D g2 = (Graphics2D) g;
                 g2.setColor(getBackground());
                 g2.fillRect(0, 0, width, height);
 
                 g2.setColor(getForeground());
                 g2.drawRect(0, 0, width, height);
                 
                 // Draw arrow and thumb backgrounds
                 g2.setStroke(arrowStroke);
                 final int thumbExtent = Math.abs(max - min) - 1;
 
                 if (orientation == VERTICAL) {
                         if (direction == LEFTRIGHT_TOPBOTTOM) {
                                 g2.setColor(getForeground());
                                 g2.fillRect(0, min - ARROW_SZ, width, ARROW_SZ - 1);
                                 paint3DRectLighting(g2, 0, min - ARROW_SZ, width, ARROW_SZ - 1);
 
                                 if (thumbColor != null) {
                                         g2.setColor(thumbColor);
                                         g2.fillRect(0, min, width, thumbExtent);
                                         paint3DRectLighting(g2, 0, min, width, thumbExtent);
                                 }
 
                                 g2.setColor(getForeground());
                                 g2.fillRect(0, max, width, ARROW_SZ - 1);
                                 paint3DRectLighting(g2, 0, max, width, ARROW_SZ - 1);
 
                                 // Draw arrows
                                 g2.setColor(arrowColor);
                                 paintArrow(g2, (width - ARROW_WIDTH) / 2.0, min - ARROW_SZ + (ARROW_SZ - ARROW_HEIGHT) / 2.0,
                                                 ARROW_WIDTH, ARROW_HEIGHT, true);
                                 paintArrow(g2, (width - ARROW_WIDTH) / 2.0, max + (ARROW_SZ - ARROW_HEIGHT) / 2.0,
                                                 ARROW_WIDTH, ARROW_HEIGHT, false);
                         } else {
                                 g2.setColor(getForeground());
                                 g2.fillRect(0, min, width, ARROW_SZ - 1);
                                 paint3DRectLighting(g2, 0, min, width, ARROW_SZ - 1);
 
                                 if (thumbColor != null) {
                                         g2.setColor(thumbColor);
                                         g2.fillRect(0, max, width, thumbExtent);
                                         paint3DRectLighting(g2, 0, max, width, thumbExtent);
                                 }
 
                                 g2.setColor(getForeground());
                                 g2.fillRect(0, max - ARROW_SZ, width, ARROW_SZ - 1);
                                 paint3DRectLighting(g2, 0, max - ARROW_SZ, width, ARROW_SZ - 1);
 
                                 // Draw arrows
                                 g2.setColor(arrowColor);
                                 paintArrow(g2, (width - ARROW_WIDTH) / 2.0, min + (ARROW_SZ - ARROW_HEIGHT) / 2.0,
                                                 ARROW_WIDTH, ARROW_HEIGHT, false);
                                 paintArrow(g2, (width - ARROW_WIDTH) / 2.0, max - ARROW_SZ + (ARROW_SZ - ARROW_HEIGHT) / 2.0,
                                                 ARROW_WIDTH, ARROW_HEIGHT, true);
                         }
                 } else { // HORIZONTAL
                         if (direction == LEFTRIGHT_TOPBOTTOM) {
                                 g2.setColor(getForeground());
                                 g2.fillRect(min - ARROW_SZ, 0, ARROW_SZ - 1, height);
                                 paint3DRectLighting(g2, min - ARROW_SZ, 0, ARROW_SZ - 1, height);
 
                                 if (thumbColor != null) {
                                         g2.setColor(thumbColor);
                                         g2.fillRect(min, 0, thumbExtent, height);
                                         paint3DRectLighting(g2, min, 0, thumbExtent, height);
                                 }
 
                                 g2.setColor(getForeground());
                                 g2.fillRect(max, 0, ARROW_SZ - 1, height);
                                 paint3DRectLighting(g2, max, 0, ARROW_SZ - 1, height);
 
                                 // Draw arrows
                                 g2.setColor(arrowColor);
                                 paintArrow(g2, min - ARROW_SZ + (ARROW_SZ - ARROW_HEIGHT) / 2.0,
                                                 (height - ARROW_WIDTH) / 2.0, ARROW_HEIGHT, ARROW_WIDTH, true);
                                 paintArrow(g2, max + (ARROW_SZ - ARROW_HEIGHT) / 2.0, (height - ARROW_WIDTH) / 2.0,
                                                 ARROW_HEIGHT, ARROW_WIDTH, false);
                         } else {
                                 g2.setColor(getForeground());
                                 g2.fillRect(min, 0, ARROW_SZ - 1, height);
                                 paint3DRectLighting(g2, min, 0, ARROW_SZ - 1, height);
 
                                 if (thumbColor != null) {
                                         g2.setColor(thumbColor);
                                         g2.fillRect(max, 0, thumbExtent, height);
                                         paint3DRectLighting(g2, max, 0, thumbExtent, height);
                                 }
 
                                 g2.setColor(getForeground());
                                 g2.fillRect(max - ARROW_SZ, 0, ARROW_SZ - 1, height);
                                 paint3DRectLighting(g2, max - ARROW_SZ, 0, ARROW_SZ - 1, height);
 
                                 // Draw arrows
                                 g2.setColor(arrowColor);
                                 paintArrow(g2, min + (ARROW_SZ - ARROW_HEIGHT) / 2.0, (height - ARROW_WIDTH) / 2.0,
                                                 ARROW_HEIGHT, ARROW_WIDTH, true);
                                 paintArrow(g2, max - ARROW_SZ + (ARROW_SZ - ARROW_HEIGHT) / 2.0,
                                                 (height - ARROW_WIDTH) / 2.0, ARROW_HEIGHT, ARROW_WIDTH, false);
                         }
                 }
                 customPaint(g2, width, height);
                 
                 g2.setColor(Color.BLACK);
                 String minS = "min", maxS = "max";
                 FontMetrics fm = getFontMetrics(g2.getFont());
                 if(numType != FLOAT)
                 {
                	 minS = Integer.toString(getLowValue());
                	 maxS = Integer.toString(getHighValue());
                 }
                 else
                 {
                	 minS = Float.toString(getLowValue()/numScale);
                	 maxS = Float.toString(getHighValue()/numScale);
                 }
            	 int sizeMax = fm.stringWidth(maxS);
                 if(min > (sizeMax+ARROW_SZ))
                     g2.drawString(minS, min-sizeMax-ARROW_SZ-5, height-2);
                 else
                     g2.drawString(minS, min+5, height-2);

                 if(max < (width - sizeMax - ARROW_SZ))
                     g2.drawString(maxS, max+ARROW_SZ+5, height-2);
                 else
                     g2.drawString(maxS, max-sizeMax-5, height-2);
         }
 
         private int pickHandle(int xOrY) {
                 int min = toScreen(getLowValue());
                 int max = toScreen(getHighValue());
                 int pick = PICK_NONE;
 
                 if (direction == LEFTRIGHT_TOPBOTTOM) {
                         if ((xOrY > (min - ARROW_SZ)) && (xOrY < min)) {
                                 pick = PICK_LEFT_OR_TOP;
                         } else if ((xOrY >= min) && (xOrY <= max)) {
                                 pick = PICK_THUMB;
                         } else if ((xOrY > max) && (xOrY < (max + ARROW_SZ))) {
                                 pick = PICK_RIGHT_OR_BOTTOM;
                         }
                 } else {
                         if ((xOrY > min) && (xOrY < (min + ARROW_SZ))) {
                                 pick = PICK_LEFT_OR_TOP;
                         } else if ((xOrY <= min) && (xOrY >= max)) {
                                 pick = PICK_THUMB;
                         } else if ((xOrY > (max - ARROW_SZ) && (xOrY < max))) {
                                 pick = PICK_RIGHT_OR_BOTTOM;
                         }
                 }
 
                 return pick;
         }
 
         public void removeChangeListener(ChangeListener cl) {
                 listeners.remove(cl);
         }
 
         public void setEmpty(boolean empty) {
                 this.empty = empty;
                 repaint();
         }
 
         public void setHighValue(int highValue) {
                 model.setExtent(highValue - model.getValue());
         }
 
         public void setLowValue(int lowValue) {
                 int e = (model.getValue() - lowValue) + model.getExtent();
                 model.setRangeProperties(lowValue, e, model.getMinimum(), model.getMaximum(), false);
                 model.setValue(lowValue);
         }
 
         public void setMaximum(int maximum) {
                 model.setMaximum(maximum);
         }
 
         public void setMinExtent(int minExtent) {
                 this.minExtent = minExtent;
         }
 
         public void setMinimum(int minimum) {
                 model.setMinimum(minimum);
         }
 
         public void setModel(BoundedRangeModel brm) {
                 model.removeChangeListener(changeListener);
                 model = brm;
                 model.addChangeListener(changeListener);
                 repaint();
         }
 
         public void setRange(int lowValue, int highValue) {
                 model.setRangeProperties(lowValue, highValue - lowValue, model.getMinimum(), model.getMaximum(),
                                 false);
         }
 
         public void setThumbColor(Color thumbColor) {
                 this.thumbColor = thumbColor;
         }
 
         protected int toLocal(int xOrY) {
                 Dimension sz = getSize();
                 int min = getMinimum();
                 double scale;
                 if (orientation == VERTICAL) {
                         scale = (sz.height - (2 * ARROW_SZ)) / (double) (getMaximum() - min);
                 } else {
                         scale = (sz.width - (2 * ARROW_SZ)) / (double) (getMaximum() - min);
                 }
 
                 if (direction == LEFTRIGHT_TOPBOTTOM) {
                         return (int) (((xOrY - ARROW_SZ) / scale) + min + 0.5);
                 } else {
                         if (orientation == VERTICAL) {
                                 return (int) ((sz.height - xOrY - ARROW_SZ) / scale + min + 0.5);
                         } else {
                                 return (int) ((sz.width - xOrY - ARROW_SZ) / scale + min + 0.5);
                         }
                 }
         }
 
         protected int toScreen(int xOrY) {
                 Dimension sz = getSize();
                 int min = getMinimum();
                 double scale;
                 if (orientation == VERTICAL) {
                         scale = (sz.height - (2 * ARROW_SZ)) / (double) (getMaximum() - min);
                 } else {
                         scale = (sz.width - (2 * ARROW_SZ)) / (double) (getMaximum() - min);
                 }
 
                 // If the direction is left/right_top/bottom then we subtract the min and multiply times
                 // scale
                 // Otherwise, we have to invert the number by subtracting the value from the height
                 if (direction == LEFTRIGHT_TOPBOTTOM) {
                         return (int) (ARROW_SZ + ((xOrY - min) * scale) + 0.5);
                 } else {
                         if (orientation == VERTICAL) {
                                 return (int) (sz.height - (xOrY - min) * scale - ARROW_SZ + 0.5);
                         } else {
                                 return (int) (sz.width - (xOrY - min) * scale - ARROW_SZ + 0.5);
                         }
                 }
         }
 
         protected double toScreenDouble(int xOrY) {
                 Dimension sz = getSize();
                 int min = getMinimum();
                 double scale;
                 if (orientation == VERTICAL) {
                         scale = (sz.height - (2 * ARROW_SZ)) / (double) (getMaximum() + 1 - min);
                 } else {
                         scale = (sz.width - (2 * ARROW_SZ)) / (double) (getMaximum() + 1 - min);
                 }
 
                 // If the direction is left/right_top/bottom then we subtract the min and multiply times
                 // scale
                 // Otherwise, we have to invert the number by subtracting the value from the height
                 if (direction == LEFTRIGHT_TOPBOTTOM) {
                         return ARROW_SZ + ((xOrY - min) * scale);
                 } else {
                         if (orientation == VERTICAL) {
                                 return sz.height - (xOrY - min) * scale - ARROW_SZ;
                         } else {
                                 return sz.width - (xOrY - min) * scale - ARROW_SZ;
                         }
                 }
         }
 
 }
