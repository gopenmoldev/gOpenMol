Attribute VB_Name = "Module1"
Public Sub tripts()
Dim app As Application
Dim dr As Drawing
Dim grs As Graphics
Dim gr As Graphic
Dim vv1 As New Vertex
Dim vv2 As New Vertex
Dim vv3 As New Vertex
Dim vvm As New Vertex
Dim vvt As New Vertex
Dim mm As Matrix
Dim xx, yy, zz, r As Double
Dim x2, y2, z2 As Double

Open "d:\tba\verfile.txt" For Input As #1
Set app = IMSIGX.Application
Set dr = app.ActiveDrawing
Set grs = dr.Graphics
Input #1, xx, yy, zz
While (xx <> -9E+99)
Set gr = grs.Add(, "TCW40SPHERE")
vv1.X = 0
vv1.Y = 0
vv1.Z = 0
gr.Vertices.AddVertex vv1
vv1.X = 0.001
gr.Vertices.AddVertex vv1
Set mm = gr.MoveRelative(xx, yy, zz)
Input #1, x2, y2, z2
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * x2, 0.1 * y2, 0.1 * z2)
Set mm = gr.MoveRelative(xx, yy, zz)
Input #1, xx, yy, zz
Wend
Close #1
End Sub
Public Sub trisphere()
Dim app As Application
Dim dr As Drawing
Dim grs As Graphics
Dim gr As Graphic
Dim vv1 As New Vertex
Dim vv2 As New Vertex
Dim vv3 As New Vertex
Dim vvm As New Vertex
Dim vvt As New Vertex
Dim mm As Matrix
Dim xx, yy, zz, r As Double
Dim x2, y2, z2 As Double

Open "d:\tba\trigfile.txt" For Input As #1
Set app = IMSIGX.Application
Set dr = app.ActiveDrawing
Set grs = dr.Graphics
For j = 1 To 4
Input #1, xx, yy, zz
vv1.X = xx
vv1.Y = yy
vv1.Z = zz
Set gr = grs.AddLineIrregularPolygon(xx, yy, zz)
Input #1, xx, yy, zz
vv2.X = xx
vv2.Y = yy
vv2.Z = zz
gr.Vertices.AddVertex vv2
Input #1, xx, yy, zz
vv3.X = xx
vv3.Y = yy
vv3.Z = zz
gr.Vertices.AddVertex vv3
Call gr.Close
Input #1, xx, yy, zz
vvm.X = xx
vvm.Y = yy
vvm.Z = zz
Input #1, r
Debug.Print r
If (r > 0) Then
    Set gr = grs.Add(, "TCW40SPHERE")
    vvt.X = 0#
    vvt.Y = 0#
    vvt.Z = 0#
    gr.Vertices.AddVertex vvt
    vvt.X = r
    gr.Vertices.AddVertex vvt
    Set mm = gr.MoveRelative(xx, yy, zz)
End If
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vv1.X, vv1.Y, vv1.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vv2.X, vv2.Y, vv2.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vv3.X, vv3.Y, vv3.Z)
Input #1, xx, yy, zz
If (r <= 0) Then
    Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Else
    frac = -r / Sqr(xx * xx + yy * yy + zz * zz)
    Set gr = grs.AddLineSingle(0, 0, 0, frac * xx, frac * yy, frac * zz)
End If
Set mm = gr.MoveRelative(vvm.X, vvm.Y, vvm.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set gr = grs.AddLineSingle(0, 0, 0, -0.1 * xx, -0.1 * yy, -0.1 * zz)
Set mm = gr.MoveRelative(vvm.X, vvm.Y, vvm.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vvm.X, vvm.Y, vvm.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vvm.X, vvm.Y, vvm.Z)
Next j
Close #1
End Sub
Public Sub trisphere2()
Dim app As Application
Dim dr As Drawing
Dim grs As Graphics
Dim gr As Graphic
Dim vv1 As New Vertex
Dim vv2 As New Vertex
Dim vv3 As New Vertex
Dim vvm As New Vertex
Dim vvt As New Vertex
Dim mm As Matrix
Dim xx, yy, zz, r As Double
Dim x2, y2, z2 As Double

Open "d:\tba\trigfile.txt" For Input As #1
Set app = IMSIGX.Application
Set dr = app.ActiveDrawing
Set grs = dr.Graphics
For j = 1 To 3
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Input #1, r
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Input #1, xx, yy, zz
Next j
Input #1, xx, yy, zz
vv1.X = xx
vv1.Y = yy
vv1.Z = zz
Set gr = grs.AddLineIrregularPolygon(xx, yy, zz)
Input #1, xx, yy, zz
vv2.X = xx
vv2.Y = yy
vv2.Z = zz
gr.Vertices.AddVertex vv2
Input #1, xx, yy, zz
vv3.X = xx
vv3.Y = yy
vv3.Z = zz
gr.Vertices.AddVertex vv3
Call gr.Close
Input #1, xx, yy, zz
vvm.X = xx
vvm.Y = yy
vvm.Z = zz
Input #1, r
Debug.Print r
If (r > 0) Then
    Set gr = grs.Add(, "TCW40SPHERE")
    vvt.X = 0#
    vvt.Y = 0#
    vvt.Z = 0#
    gr.Vertices.AddVertex vvt
    vvt.X = r
    gr.Vertices.AddVertex vvt
    Set mm = gr.MoveRelative(xx, yy, zz)
End If
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vv1.X, vv1.Y, vv1.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vv2.X, vv2.Y, vv2.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vv3.X, vv3.Y, vv3.Z)
Input #1, xx, yy, zz
If (r <= 0) Then
    Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Else
    frac = -r / Sqr(xx * xx + yy * yy + zz * zz)
    Set gr = grs.AddLineSingle(0, 0, 0, frac * xx, frac * yy, frac * zz)
End If
Set mm = gr.MoveRelative(vvm.X, vvm.Y, vvm.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set gr = grs.AddLineSingle(0, 0, 0, -0.1 * xx, -0.1 * yy, -0.1 * zz)
Set mm = gr.MoveRelative(vvm.X, vvm.Y, vvm.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vvm.X, vvm.Y, vvm.Z)
Input #1, xx, yy, zz
Set gr = grs.AddLineSingle(0, 0, 0, 0.1 * xx, 0.1 * yy, 0.1 * zz)
Set mm = gr.MoveRelative(vvm.X, vvm.Y, vvm.Z)
Close #1
End Sub

Public Sub drawstuff()
Call tripts
Call trisphere2
End Sub
