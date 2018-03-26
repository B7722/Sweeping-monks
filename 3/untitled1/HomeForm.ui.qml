import QtQuick 2.7
import QtQuick.Controls 2.0

Page {
    width: 600
    height: 440
    title: qsTr("Home")

    Image {
        id: image
        x: 0
        y: 0
        width: 600
        height: 440
        fillMode: Image.Stretch
        source: "images/timg.jpg"
    }
}
