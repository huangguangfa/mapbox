# mapbox 学习旅途
mapbox学习

在学习mapbox的时候先了解什么是GIS，可以看下这个文章`https://zhuanlan.zhihu.com/p/479427335?utm_medium=social&utm_oi=53772994740224`

## 文档
[mapbox官网](https://www.mapbox.com)

[mapbox文档(JS)](https://docs.mapbox.com/mapbox-gl-js)

[GeoJson 格式规范](https://www.rfc-editor.org/rfc/rfc7946.html)

[Turf(地理空间分析库)](http://turfjs.org)

## 本仓库结构介绍

```js
- mapbox-gl-js // mapbox开源渲染引擎
- example // 案例

```

> 案例的类型主要是以layers文档类型来，source来源可以参考`GeoJson的规范`, 另外阅读的方式提高自己想学习的模块来进行学习

## 开始
```js
cd mapbox-gl-js   // 进入maobox 源码进行启动
nvm use 16        // 切换node 16版本
yarn install      // 安装依赖
yarn start        // 启动项目
```
正常看到maobox-gl-js里面有dist产物就算可以了，我们只要产物就行