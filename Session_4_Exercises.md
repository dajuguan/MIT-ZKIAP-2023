# Set up a starter project
按照教程走即可，见`my-app`文件夹

# [Boilerplate](https://github.com/semaphore-protocol/boilerplate)
这部分其实官方[Demo](demo.semaphore.appliedzkp.org/)最后一步的**Proofs**无法运行，这是由于跨域了，可以采用[CORS代理](https://stackoverflow.com/questions/43262121/trying-to-use-fetch-and-pass-in-mode-no-cors)的方式来增加相应的Header。
所以需要修改`boilerplate/apps/contracts/scripts/download-snark-artifacts.ts`文件中的`url`为:
``` js
// const url = `http://www.trusted-setup-pse.org/semaphore/${process.env.TREE_DEPTH}`
const url = `https://cors-anywhere.herokuapp.com/https://www.trusted-setup-pse.org/semaphore/${process.env.TREE_DEPTH}`
```
不想自己编译，就体验一下的话，可以采用浏览器override的方式，将`/demo.semaphore.appliedzkp.org/_next/static/chunks/120-8439565389d7b71e.js`中相应语句修改为
```
return snarkArtifacts || (snarkArtifacts = {
    wasmFilePath: "https://cors-anywhere.herokuapp.com/https://www.trusted-setup-pse.org/semaphore/".concat(merkleProof.siblings.length, "/semaphore.wasm"),
    zkeyFilePath: "https://cors-anywhere.herokuapp.com/https://www.trusted-setup-pse.org/semaphore/".concat(merkleProof.siblings.length, "/semaphore.zkey")
}),
```