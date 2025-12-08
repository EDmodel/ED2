# How to contribute to the ED2 model

Thank you for your interest in contributing to the ED2 model development. ED2 development takes a community-based approach, so anyone working with the model is very welcome to contribute. We only ask you to follow some contributing guide lines to ensure code quality and reproducibility.

## Before you start

Please make sure to read, understand and accept the following documents:

- **ED2 Licence**. The ED2 model is licenced under the _Creative Commons Attribution 4.0 International_ (CC-BY-4.0). You can find the complete licence in [this link](https://github.com/EDmodel/ED2/blob/master/LICENSE).

- **ED2 Code of Conduct**. To ensure everyone willing to collaborate with the ED2 model development can do so in a positive environment, we ask everyone to follow our [code of conduct](https://github.com/EDmodel/ED2/blob/master/CODE_OF_CONDUCT.md).


## Steps for contributing with code development

To make sure your contributions are the most effective, please consider the following steps:

1. **Communicate your intent to the community ahead of time.** If you would like to propose a new feature, or if you think you found a bug in the code, the best way to start is by checking previous [issues](https://github.com/EDmodel/ED2/issues) and [discussions](https://github.com/EDmodel/ED2/discussions). If none of them seem to address your proposed change, then [open a new issue](https://github.com/EDmodel/ED2/issues) that describes your plans. This allows the community to engage with your plans at an early stage and may avoid redundant development.
2. **Create a branch from the main branch**. The official repository of ED2 is [https://github.com/EDmodel/ED2](https://github.com/EDmodel/ED2). We **strongly** recommend using this repository as your starting point and create a fork from this repository. This will help code reviewers to efficiently give feedback to your contributions and maintainers to merge your contributions.
3. **Less is more**. Smaller and focussed contributions to the model are much preferred, because they are easier to understand and to handle. That said, we understand that sometimes the model development may require more substantial changes. If this is the case, one helpful approach is to make multiple [git branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches), each containing one small and focussed component of a broader development.
4. **Documentation is fundamental**. Documentation comes in multiple places, and they are all important. First, make sure that your changes in the code come with comments in the code itself. When in doubt, err on the side of explaining more, it will help others understand the changes, and why they were needed. Second, make sure to provide as much context as possible when submitting a pull request (see step 7). If you implement new features, it is important to update the [ED2 Wiki](https://github.com/EDmodel/ED2/wiki) to explain the features.
5. **Style is fundamental too**. Keeping with the best coding practices and style makes the code more readable and robust, and will help reviewers and maintainers to merge your contributions. Please refer to our [code organisation and design](https://github.com/EDmodel/ED2/wiki/Code-organization-and-design-philosophy) best practices.
6. **Test, test, test**. Before submitting your contribution, make sure you test your changes. One very useful step is to compile the code and run it with comprehensive and strict debugging flags, which may help you find some code bugs at an early development stage. Please refer to our [compilation instructions](https://github.com/EDmodel/ED2/wiki/Compiler-instructions-%28aka-the-include.mk-files%29) for more details on how to build ED2 with strict compilation rules.
7. **Submit a pull request**. Once you are confident that your changes are ready, create a [pull request](https://github.com/EDmodel/ED2/pulls). Make sure to document your changes in detail, and engage with the community as you receive feedback. At this point, it is very common that people will provide you suggestions on additional modifications. Please follow up with them, and in most cases, it is possible to submit additional changes to the same pull request. Once the changes are approved, your contributions will be merged to the ED2 main branch.

## Additional resources

- [GitHub Docs](https://docs.github.com/en). An excellent starting point for best practices when using git and GitHub for ED2 or any other project. Particularly relevant documents include (1) [how to use forks](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks), (2) [how to use branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches) and (3) [how to make pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).
- [ED2 Wiki](https://github.com/EDmodel/ED2/wiki). This is our live documentation for ED2. Like the model itself, the documentation is also community-based, so feel free to contribute to it if any information is missing or needs update.

